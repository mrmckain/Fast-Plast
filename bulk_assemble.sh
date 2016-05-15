#!/bin/bash

cd `pwd`
CURRDIR=`pwd`
DATA=$2
while read GENE
do
        echo $GENE
        NAMESP=$(echo $GENE | cut -d " " -f 2)
	TKNUM=$(echo $GENE | cut -d " " -f 1)
        echo $TKNUM
	NEWDIR=$NAMESP
        NEWDIR+="_"
        NEWDIR+=$TKNUM
        REAL=$(echo $NEWDIR | cut -d " " -f 2)
	echo $REAL
	mkdir $REAL
        cd $REAL
        cp ../plastome-assembly_condor_submit .
       	FILELIST=("$DATA/$TKNUM"*fastq)
	FILE1=${FILELIST[0]}
	echo $FILE1
	FILE2=${FILELIST[1]}
        echo $FILE2
	
	NEWNAME="$FILE1 $FILE2 $NAMESP""_""$TKNUM"
	echo $NEWNAME
	cat  plastome-assembly_condor_submit | sed -e "s|REPLACE|$NEWNAME|g" > plastome-assembly_condor_submit_temp
	chmod +x plastome-assembly_condor_submit_temp 
        condor_submit plastome-assembly_condor_submit_temp
	cd $CURRDIR
done < $1
