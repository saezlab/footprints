#!/bin/bash

cd $(dirname $0)

rsync -auvr --include='*.pdf' --include='*/' --exclude='*' \
    ebi:/nfs/research2/saezrodriguez/mike/speed2/* . 
