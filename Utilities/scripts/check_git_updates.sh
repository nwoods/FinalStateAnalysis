#! /bin/bash

echo $PWD
#check that on-line repository uwcms is not ahead of us
remoteName=$(git remote -v | grep uwcms | awk '{print $1}' | head -n 1) #finds which name we gave to uwcms repository
lastCommit=$(git ls-remote $remoteName | grep HEAD | awk '{print $1}') #gets the latest commit hash in the HEAD of UWCMS

if [ -z $(git log | grep $lastCommit | awk '{print $2}') ]; then #looks if we have that commit in our git log (the print $2 is just to get rid of the word 'commit' that gives bash some problem
    echo "UWCMS master repository is ahead of the used one! Should I pull the changes? [y/n]"
    read choice
    if [ $choice = 'y' ]; then
	echo "pulling from $remoteName master"
	git pull $remoteName master
	echo "remember to recopile! [any key]"
	read ok
    else
	echo "ok, but remember to do it asap! [any key]"
	read ok
    fi
else
    echo "This area is up do date!"
fi

#ok, now that we know we are up to date check if WE have someting to push in uwcms
commitsAhead=$(git log | grep commit | grep -n $lastCommit | awk -F: '{print $1}')

if [ -n $commitsAhead ]; then
    if [ $commitsAhead != '1' ]; then
	echo "Your repository is ahead of $commitsAhead, consider sending a Pull Request! [any key]"
	read ok
    fi
fi
