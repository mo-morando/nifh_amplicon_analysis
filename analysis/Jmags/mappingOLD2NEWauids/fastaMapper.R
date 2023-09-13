#!/usr/bin/env Rscript

usageStr="Usage:
    fastaMapper fasta1 fasta2

Write to fasta_id_map.tsv the IDs of sequences that exactly match between two FASTA files in the
following format:
    [ID in fasta1]<tab>[ID in fasta2]
where ID is everything after the '>' up to the first space.  Only IDS for sequences that appear in
both files identically (ignoring case) will be output.

FASTA sequences may be split over multiple lines.  Blank lines can appear between records.

This script is simple.  It assumes that each sequence in fasta1 is unique and has a unique
definition line, and similar for fasta2.  Results are undefined if a sequence is repeated within
the FASTA (with distinct or repeated definition lines).
"

args <- commandArgs(T)
if (! (file.exists(args[1]) && file.exists(args[2]))) {
    cat(usageStr)
    stop("Need two fasta files.")
}

cat("Reading in ShortRead lib...")
suppressMessages(library(ShortRead))
cat("done.  Hold on while I make fasta_id_map.tsv...")

fas1 <- readFasta(args[1])
fas2 <- readFasta(args[2])
shared <- intersect(sread(fas1), sread(fas2))
idx1 <- match(shared, sread(fas1));  stopifnot(!is.na(idx1))
idx2 <- match(shared, sread(fas2));  stopifnot(!is.na(idx2))
dat <- data.frame(ids1 = sub(' .*', '', id(fas1)[idx1]),
                  ids2 = sub(' .*', '', id(fas2)[idx2]))
write.table(dat, "fasta_id_map.tsv", sep="\t", quote=F, row.names=F)
cat("done.\n")
quit(save='no')
