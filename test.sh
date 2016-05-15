cd `pwd`
CURRDIR=`pwd`
while read GENE
do
        #echo $GENE
        NAMESP=$(echo $GENE | cut -f 2)
        TKNUM=$(echo $GENE | cut -f 1)
        NEWDIR=$NAMESP
	NEWDIR+="_"
	NEWDIR+=$TKNUM
        REAL=$(echo $NEWDIR | cut -d " " -f 2)
	echo $REAL
done < $1
