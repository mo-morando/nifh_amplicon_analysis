#!/usr/bin/bash

##
## Map sequence identifiers between two FASTA files.
##

usageStr="Usage:
    fastaMapper fasta1 fasta2 [options]

Write to stdout the IDs of sequences that exactly match between two FASTA files in the following
format:
    [ID in fasta1]<tab>[ID in fasta2]
where ID is everything after the '>' up to the first space.  Only IDS for sequences that appear in
both files identically (ignoring case) will be output.

Options:
    -seqs      Include a third column in the output for the sequence.
    -help      Help!
    -h

FASTA sequences may be split over multiple lines.  Blank lines can appear between records.

This script is simple.  It assumes that each sequence in fasta1 is unique and has a unique
definition line, and similar for fasta2.  Results are undefined if a sequence is repeated within
the FASTA (with distinct or repeated definition lines).
"
if [ ! -z `echo "$@" | egrep "(-h|-help)"` ] ; then
    echo "$usageStr"
    exit 0
fi

fas1=$1
fas2=$2
if [ ! -f "$fas1" ] || [ ! -f "$fas2" ] ; then
    echo "ERROR: Both fasta files must exist."
    exit -1
fi

## Make sequences appear on just one line ("linearize").  Code from here:
##  https://www.biostars.org/p/17680/#226890
## Also removes empty lines.
LinearizeFasta()
{
    awk '{ if (NR==1) {
               print $0
           } else {
               if ($0 ~ /^>/) { print "\n"$0 } else { printf $0 }
           }
         }' "$1" > "$2"
}

## Make a fasta into a table with col 1 the ID and col 2 the sequence.
## Sort the table alphabetically by the sequence column. (Needed for 'join'.)
Fasta2Table()
{
    fas=$1
    outfile=$2
    LinearizeFasta "$fas" fas.tmp
    cat fas.tmp | grep '^>' | sed 's/^>\([^ ]*\).*$/\1/' > ids.tmp
    cat fas.tmp | grep -v '^>' > seqs.tmp
    paste ids.tmp seqs.tmp | sort -k 2b,2 > "$outfile"
    rm fas.tmp seqs.tmp ids.tmp
}

Fasta2Table "$fas1" fas1.tab
Fasta2Table "$fas2" fas2.tab
if [ -z "`echo $@ | grep "-seqs"`" ] ; then
    join -t $'\t' -1 2 -2 2 fas1.tab fas2.tab | cut -f2,3
else
    join -t $'\t' -1 2 -2 2 fas1.tab fas2.tab | awk '{print $2"\t"$3"\t"$1}'
fi
rm fas1.tab fas2.tab
exit 0
