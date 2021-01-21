#!/usr/bin/bash
#SBATCH -p short
module load samtools/1.11
module load bwa/0.7.17
if [ -f config.txt ]; then
	source config.txt
fi
mkdir -p $GENOMEFOLDER
pushd $GENOMEFOLDER
FASTAFILE=candida_lusitaniae_ATCC42720_w_CBS_6936_MT.fasta
GFF=candida_lusitaniae_1.sorted.gff3
## THIS IS FUNGIDB DOWNLOAD PART
echo "working off $FASTAFILE - check if these don't match may need to update config/init script"

if [[ ! -f $FASTAFILE.fai || $FASTAFILE -nt $FASTAFILE.fai ]]; then
	samtools faidx $FASTAFILE
fi
if [[ ! -f $FASTAFILE.bwt || $FASTAFILE -nt $FASTAFILE.bwt ]]; then
	bwa index $FASTAFILE
fi

DICT=$(basename $FASTAFILE .fasta)".dict"

if [[ ! -f $DICT || $FASTAFILE -nt $DICT ]]; then
	rm -f $DICT
	samtools dict $FASTAFILE > $DICT
	ln -s $DICT $FASTAFILE.dict 
fi

popd
