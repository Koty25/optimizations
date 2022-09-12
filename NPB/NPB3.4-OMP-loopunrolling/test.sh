#!/bin/bash

nlines=10
nlines="$(($nlines-1))"
for line in  $(seq 1 $nlines); do
echo $line
done
