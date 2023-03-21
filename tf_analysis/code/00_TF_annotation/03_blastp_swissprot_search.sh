# Using the previous blastp output from swissprot, we blast swissprot against our set of sequences of interest, using a stringent pvalue.
# If this output is not detected, then we generate a blast from the sequence of interest against uniprot.
# From this, we grep the words "Transcription factor", "transcription factor", to retrieve the matches that are transcription factors.
# We retrieve the columns with the same matches. This is what we take as reciprocal best hits.

# Parameters: pep.fa, the swissprot database
# output: the tsvs of blast, and the 2-column table of protein sequence ids <--> names of hits