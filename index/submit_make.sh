#!/bin/bash

test $1 || exit "need to supply target and/or options"

bsub -M 1024 -R "rusage[mem=1024]" "cd $(dirname $0) && make $*"


