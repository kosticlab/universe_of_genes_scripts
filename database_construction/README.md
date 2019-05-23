## Directory overview

The folder contains the scripts we used to generate our database, which is located at microbial-genes.bio.

## Directory contents

#### format_db.py

This is the script we used to generate the database file from our various data files (gene taxonomic annotations, gene catalogs, functional annotations). Note that in order to compute the number of genes in each cluster (i.e. identify singletons vs. non-singletons), it was necessary to track which genes clustered with what across the oral and gut catalogs. This is because each catalog was made individually before being merged together.

#### format_synthetic_db.py

Same as above, but we used it to test our database generation at small scale.

#### website_plots.R

An R script containing code to generate the gene length plot found on our website.
