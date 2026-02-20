A python script that (1) filters a predicted-terminator CSV by a settable --threshold on probability_mean (default 0.3), (2) BLASTs the FASTA queries against a local BLAST DB, (3) computes percent overlap between BLAST-hit genomic intervals and predicted terminator intervals in the same genome/contig, and (4) computes a random-region control plus an optional plot. 
The BLAST parsing uses tabular -outfmt 6 fields including qseqid, sseqid, qstart, qend, sstart, send, evalue, bitscore as documented for BLAST+ tabular output.
​

Usage
# Basic
python terminator_overlap.py \
  --pred-csv predicted_terminators.csv \
  --query-fasta queries.fa \
  --blast-db /path/to/local_blastdb/prefix \
  --out-prefix results/out \
  --threshold 0.3

# With plot + more random controls
python terminator_overlap.py \
  --pred-csv predicted_terminators.csv \
  --query-fasta queries.fa \
  --blast-db /path/to/local_blastdb/prefix \
  --plot \
  --random-n 1000 \
  --out-prefix results/out


dependencies:
pandas
biopython
intervaltree
blast
matplotlib