#!/bin/bash

cd $(dirname $0)

rsync -auvr --include='*.pdf' --include '[zd]scores.RData' --exclude 'data/*' \
    --include '*.RData' --include='*/' --exclude='*' \
    ebi:/nfs/research2/saezrodriguez/mike/speed2/* . 
