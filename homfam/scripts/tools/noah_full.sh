#!/bin/bash
echo $1
echo $2
./ramsaics -i "$1" -o "$2" --arc-align --phylo
