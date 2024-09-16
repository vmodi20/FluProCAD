#!/bin/bash

mkdir slow-TI

for i in $(seq 0 1 20); do
	mkdir slow-TI/$i
	for j in *fe.mdp setup-feMD.sh ; do
		sed "s/@lambda@/$i/" $j >& slow-TI/$i/$j
	done
done
