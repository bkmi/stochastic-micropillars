#!/bin/bash

gss=$1
gww=$2
gws=$3
gsw=$4
rawargs='-sponE -makeavgs 20 -t_stop 10000 -rmlog -corrprintlen 40 -fourierprintlen 40 -fA_ss 0.0 -fA_ww 0.0 -tau_fb 0'
mem=6



#for i in `seq 0 1 2`; do
#	spath=' -spath changetau_'$i'_FB005/ '
#	args=$rawargs$spath' -tau_fb '$i
#	echo $args
#	qsub -mem $mem -m n -o /dev/null -args "$args" muPillar_FB
	#qsub -mem $mem -m n  -args "$args" muPillar_FB
#	done


