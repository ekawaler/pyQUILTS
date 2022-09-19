#!/bin/bash

<<COMMENT
    For each .vcf file in the directory:
        Turn it into a .vcf.gz file
    For each .vcf.gz file in the directory:
        Reindex, or it throws a tantrum
        Check whether it's a multi-sample or single-sample file
            If single: add format column & sample column with name of file as sample name, along with GTs
        Check chromosome name format
            If incorrect: Rename chromosomes according to format file [error if it's a format or chromosome not found]
            At the moment, can handle format "chr#" or just "#"...explain this whole thing better
        Unroll multiple alleles
    Merge
    TO DO: stick quality filter in there somewhere. Use QUAL or GQ?
        Will look something like this: bcftools view -i '%QUAL>=15' gtypes.vcf.gz
COMMENT

FILES=( "$1"/*.vcf )
if [ -f "$FILES" ]
then
    for f in "${FILES[@]}"
    do
        bgzip -c $f > "${f}".gz
    done
else
    :
fi

FILES=( "$1"/*.vcf.gz )
if [ -f "$FILES" ]
then
    for f in "${FILES[@]}"
    do
        tabix -p vcf $f
        bcftools norm -m - $f -o tmp.vcf # unrolls multi-allele SNPs
        bcftools annotate --rename-chrs "chrmap.txt" tmp.vcf -o tmp2.vcf
        bgzip -c tmp2.vcf > $f
        tabix -p vcf -f $f
        rm tmp.vcf
        rm tmp2.vcf
        samps=$(bcftools query -l $f | wc -l)
        if [ $samps == "0" ]
        then # Add format/sample columns
            filname=$(basename -- "$f")
            sampname=${filname%".vcf.gz"}
            bcftools view -h $f | grep '^##' > header.hdr # to get what we need from the header
            echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> header.hdr
            echo $'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'$sampname >> header.hdr
            bcftools view -H $f | sed '/^#/!s/$/\tGT\t1\/\./' > tmp.vcf
            cat header.hdr tmp.vcf > tmp2.vcf
            bgzip -c tmp2.vcf > $f
            tabix -p vcf -f $f
            rm tmp.vcf
            rm tmp2.vcf
        fi
    done
    if [ ${#FILES[@]} -gt 1 ]
    then
        bcftools merge -0 -m none -o $1/filtered.vcf "${FILES[@]}"
    else
        bcftools view $f > $1/filtered.vcf
    fi
else
    exit 1
fi
