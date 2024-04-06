#!/bin/bash
set -e
set -u
set -o pipefail

sample_name=$(basename -s ".txt" "$1")
awk '{{gsub(/+/, "F");gsub(/-/,"R"); print $1"."$2, $1, $2, $3, $6=$4+$5, $7=100*$4/$6, $8=100*$5/$6}}' ${sample_name}.txt \
| awk 'BEGIN{print "chrBase	chr	base	strand	coverage	freqC	freqT"}1' \
| awk 'BEGIN{OFS ="\t"}{print $1, $2, $3, $4, $5, $6, $7}' > 010c.methylkit/${sample_name}-MK.txt