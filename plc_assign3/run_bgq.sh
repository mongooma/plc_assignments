#!/bin/sh
#salloc --nodes 256 --time 60 --partition large # not necessary when not testing
mpixlc -O5 main.c functions.c -o main.xl