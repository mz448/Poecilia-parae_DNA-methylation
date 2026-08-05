#!/bin/bash

# concatenate all features in a single file
cat *.feature.bed > allFeatures.feature.bed

# Same using a for loop
# for FEATFILE in *.feature.bed; do
#	cat $"{FEATFILE}" > allFeatures.feature.bed
#	done
