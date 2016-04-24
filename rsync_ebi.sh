#!/bin/bash

cd $(dirname $0)

rsync -auvr --include='*.pdf' --include '[zd]scores.RData' --exclude 'data/*' \
    --include '*.RData' --include='*/' --exclude='*' \
    yoda:/research/performance/saezrodriguez/mike/speed2/* . 

#rsync -avur figure/* yoda:/nfs/research2/saezrodriguez/mike/speed2/figure/
