#!/bin/bash

for i in `find DNA_alignments -type f -print0 | gshuf | head -n 300`
do
	echo $i
done
