#/!bin/bash

SMRTDIR=$1

# link filtered_subreads.fasta
fa=$SMRTDIR/data/filtered_subreads.fasta

if [ ! -e $fa ] 
then
    echo "$fa does not exist! ABORT!"
    exit
fi

ln -s $fa .
cp ../primers.fa .

echo """
#!/bin/bash

. ~/.VENV/bin/activate
PacBioBarcodeIDCCS.py filtered_subreads.fasta ../primers.fa output || exit $?;
barcode_trimmer.py -i filtered_subreads.fasta -d output --output-anyway -o filtered_subreads.trim_both_anyway.fa > filtered_subreads.trim_both_anyway.fa.log || exit $?;
seqclean filtered_subreads.trim_both_anyway.fa -c 8 || exit $?;
parse_seqclean_cln.py filtered_subreads.trim_both_anyway.fa.report filtered_subreads.trim_both_anyway.fa.cln || exit $?;
"""
