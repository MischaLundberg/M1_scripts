#!/bin/sh

if [ $# -ne 3 ]
then
    echo "usage: $0 <GWAS INPUT FILE> <ALIAS> <SAMPLE SIZE>"
    exit 65
fi


### Setup preferences ###
#python=/mnt/bin/anaconda2/bin/python2.7
python=python
ldsc_dir=/mnt/data/ldsc/

gwas=$1
alias=$2
N=$3

Rscript bin/check_integrity.R $gwas ${gwas}.checked 1200 700 2 >${gwas}.checked.log 2>${gwas}.checked.err || { echo "Integrity check failed"; exit 1; }

$python $ldsc_dir/munge_sumstats.py --sumstats ${gwas}.checked --N $N --out ${gwas}.checked --merge-alleles $ldsc_dir/w_hm3.snplist >${gwas}.munge.log 2>${gwas}.munge.err || { echo "Formatting failed"; rm $gwas*; exit 1; }

#echo Chr ProbeID GeneticDistance ProbeBp Gene Orientation PathOfEsd >${gwas}.flist
#echo 1 $alias 1 1 $alias + ${gwas}.checked >>${gwas}.flist

#./bin/smr_3 --eqtl-flist ${gwas}.flist --make-besd --out ${gwas} --peqtl-trans 1 --peqtl-other 1 || { echo "Indexing failed"; exit 1; }

#if [ ! -f ${gwas}.besd ]; then
#	echo "Problem indexing"
#	exit 1
#fi

bgzip ${gwas}.checked 2>${gwas}.bgzip.err || { echo "Compression failed"; rm $gwas*; exit 1; }
tabix -s 1 -b 3 -e 3 -S 1 ${gwas}.checked.gz 2>${gwas}.tabix.err || { rm $gwas*; echo "Indexing failed"; exit 1; }

rm $gwas

