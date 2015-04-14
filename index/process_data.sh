#!/bin/bash

(for DIR in $(ls | egrep -o "GSE[0-9]+"); do
    echo $DIR
    cd $DIR
    for FILE in $(ls *.tar); do 
        tar -xvf $FILE
    done
    gunzip *.gz
    rm *.tar 
    cd ..
done )>ext.out 2>ext.err

