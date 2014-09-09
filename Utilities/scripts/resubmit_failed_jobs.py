#!/usr/bin/env python
"""
Allows the resubmission of failed Condor jobs provided they have a rescue dag
file.

Author: D. Austin Belknap, UW-Madison
"""

import os
import re
import sys
import glob
import argparse


def submit_jobid(jobid, dryrun=False):
    """
    Scan through the samples of a given job id, and check the dag status file
    for failed jobs. If any, submit the rescue dag files to farmoutAnalysisJobs
    """
    user = os.environ['USER']
    scratch = '/nfs_scratch'

    path = os.path.join(scratch, user, jobid)

    samples = glob.glob('%s/*' % path)

    for s in samples:
        # FSA ntuples and PAT tuples use a different naming convention for the
        # status dag files. Try both.
        status_dag1 = '%s/dags/dag.status' % s
        status_dag2 = '%s/dags/dag.dag.status' % s

        # look for failed jobs
        errors = []
        try:
            with open(status_dag1, 'r') as dagfile:
                errors = [re.search('STATUS_ERROR', line) for line in dagfile]
            with open(status_dag1, 'r') as dagfile:
                submitted = [re.search('STATUS_SUBMITTED', line) for line in dagfile]
        except IOError:
            try:
                with open(status_dag2, 'r') as dagfile:
                    errors = [re.search('STATUS_ERROR', line) for line in dagfile]
                with open(status_dag2, 'r') as dagfile:
                    submitted = [re.search('STATUS_SUBMITTED', line) for line in dagfile]
            except IOError:
                print "Skipping: %s" % s
                continue

        # Do not try to resubmit jobs if jobs are still running
        if any(submitted):
            print "'%s' is still running -- wait for jobs to complete" % s
            continue

        # if there are any errors, submit the rescue dag files
        if any(errors):
            if dryrun:
                print "Resubmit: %s" % s
            else:
                rescue_dags = glob.glob('%s/dags/*dag.rescue[0-9][0-9][0-9]' % s)
                cmd = 'farmoutAnalysisJobs --rescue-dag-file=%s' % max(rescue_dags)
                os.system(cmd)


def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Resubmit failed Condor jobs')

    parser.add_argument('jobid', type=str, help='Provide the FSA job ID the'
                                                'original jobs were run with')

    parser.add_argument('--dry-run', dest='dryrun', action='store_true',
                        help='Show samples to submit without submitting them')

    args = parser.parse_args(argv)

    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    submit_jobid(args.jobid, dryrun=args.dryrun)

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
