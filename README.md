# protein_properties

'data' directory requires 'ecod.develop279.domains.txt' and 'ecodf.hmm' (ECOD v279).

'genome_data' directory requires HMMsearch outputs of ECOD v279 domains against genomes from archaea/bacteria (GTDB v207) and eukaryotes (EukprotV3) with corresponding file names.
Taxonomy files are uploaded to Github, but not the HMMsearch outputs (size limit).

- main.py for calculating Phyletic Distribution Scores of ECOD domains in aforementioned genomes
- aerobe_annotation.py for analyzing X-group distributions in bacteria genomes with oxygen physiology annotations
