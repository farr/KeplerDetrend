#!/bin/bash

set -e 

for q in {01..17}; do
    mkdir -p HLSP/Q${q}
    mv *-q${q}_*.fits HLSP/Q${q}/
    echo "Moved files for Q${q}"
done