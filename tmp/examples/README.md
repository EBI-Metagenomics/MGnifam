TEST EXAMPLES

This folder contains manually retrieved SEED alignments and HMMs, retrieved by the first 10 clusters of batch number 100, running:

`head -n 10 build_list_0100 | awk 'print {$1}' > batch_list_0100`

First, `mgseed.pl` has been run for every of these clusters, obtaining SEED alignments

Then, `pfbuild -withpfmake -db uniprot` has been run inside each cluster folder, creating the HMM from given SEED alignment and comparing it with UniPRot .fasta sequences dataset.
