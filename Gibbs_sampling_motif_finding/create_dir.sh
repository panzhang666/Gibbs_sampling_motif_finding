#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
i=0
j=1
while [ $j -ge 1 ]
do 
    echo "Welcome $i times"
    i=$(bc <<< "$i+0.1")
    echo "$i"
    j=$(bc <<< "$i <=2")
    echo "$j"
done
