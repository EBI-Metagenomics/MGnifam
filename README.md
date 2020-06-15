# MGnifam

MGnifam dataset contains Domains of Unknown Function (DUFs) not contained in
Pfam dataset. This projects aims at the development of an automated pipeline
for the creation of a MGnifam release.

## The pipeline

The pipeline is exploited in detail at `schema/schema.graphml` (the schema has
been created using yEd diagramming library, freely available either to
download and either as live verision at https://www.yworks.com/yed-live/)

### 1. Clustering and batch creation

Input of the pipeline is the union of two sequences sets: the previous version
of MGnify (if any) and UniProt, then clusterized using LinCLust algorithm.

However, clusters are considered as input of the pipeline, since LinClust
is assumed to be run previously.

Clusters are then grouped together in batches, each different batch undergoes
the pipeline separately (until a few step before the release, when all the
batches together must be checked for overlaps).

*TODO* Integrate LinClust step inside the pipeline.

### 2. Seed creation

The first step in processing a single batch of clusters is the creation of the
SEED alignment of all the sequences contained in the cluster itself. Actually,
the script `mgseed.pl -cluster <cluster name>` handles this step, creating a
new `<cluster name>/` directory for each cluster.

Note that `mgseed.pl` either searches for intrinsically disordered regions
which, if found, makes the cluster to be discarder, i.e. put into `BIAS/`
directory.

*TODO* Implement `mgseed.pl` as a python module, integrate it with bsub.

### 3. HMM creation

Afterwards, from each seed alignment previously created, an HMM is developed
through `cd <cluster name>/ && pfbuild -withpfmake -db uniprot` command.

Once the HMM is in place, it must be run against Uniprot to check if matches
any sequence in it, which would cause the current cluster to be discraded,
i.e. put into `./Uniprot` since, it is a possible PFam and not a possible
MGnifam anymore.

### 4. Build initialization

In this step, all clusters remaining in initial batches, i.e. the ones not
moved in `BIAS/` or `Uniprot/` directory are put all together in a single
directory, here called `Build/`.

For every cluster in `Build/` directory, seed alignment must either be trimmed
to N- and C- terminal and then made non redundant, in order to avoid keeping
useless information stored on disk.

As in step 3., `cd <cluster name>/ && pfbuild -withpfmake -db uniprot` is run
after seed trimming.

Since, with very high probability, seed alignment changed, then either HMM
has to be built again as in step 3., by running
`cd <cluster name>/ && pfbuild -withpfmake -db uniprot`.

Either the instersection of clusters HMMs with Uniprot dataset is not ensured
to be null anymore, hence `check_uniprot.pl` script has to be run again,
which ensures null intersection again and discards possible PFam ones.

Note that mean batch size should be in the order of 1/10 with respect to
initial batch size. However, there will probably be more than 10 clusters,
hence the `Build/` directory will probably be (much) greater than a randomly
picked batch.

*TODO* Automatic seed trimming must be integrated.

*TODO* Non redundancy must be integrated.

### 5. Build against MGnifam

In this step, remaining HMMs (i.e. the ones which do not have any overlap with Uniprot and no low disordered regions) are run against MGnify by executing
`cd <cluster name>/ && mfbuild`

Note that `mfbuild` is just an alias for `pfbuild -db mgnify -withpfmake -makeEval 0.01`

### 6. Annotation

`DESC` file must be write out in order to define some description
attributes, such as annotiation author. To do this first a new MGnifam accession
number must be refrieved, by running `nextMGDUF.pl`.

Once the MGnifam accession number has been retrieved, one can edit `DESC` file manually by setting:
1. ID (identifier) line as `ID <new MGDUF>`
2. AU (author) line as `AU <author name>`
3. DE (description) line as `DE Protein of unknown function (<new MGDUF>)`

Otherwise, one can simply run `annotateMGDUF.pl -directory . -author '<author name>'`
insde the `Build/` directory

### 7. Databade updating

Clusters which have not been excluded from `Build/` directory, must then be
moved to a new `checked_in/` directory (i.e. possible Pfam must be discarded),
then `mfnew.pl` script must be run in each of the cluster's directories, by running
`cd checked_in/<cluster name>/ && mfnew.pl`

`mfnew.pl` checks in clusters into MySQL datbase and assigns them a new MGnifam accession number.

### 8. Generating flatfiles

Flatfiles generation is the last step before producing a new MGnifam release.
However, this step undergo a couple of checks before generating those files.

First control step checks if families in MGnifam are related to each other
and/or to families in Pfam, by generating a consensus sequence for each MGnifam
entry (using `hmmemit`) and running this against the MGnifam and Pfam models
already in place. Moreover, consensus sequences are created for Pfam entries
either and run against MGnifam entries.

The script used to generate consensus sequences and make comparisons
automatically is
`mgnifam_pfam_comparison.pl -pfam_dir <path/to/Pfm-A.hmm> -pfam_config <path/to/pfam/config> -mgnifam_config <path/to/mgnifam/config>`

Second control steps checks for midpoint overlaps using Jaccard index and
jaccard containment values for each pair of families in MGnifam. The script
used to accomplish this evaluation is `jaccard.pl`

*TODO* `jaccard.pl` does not scale well and must be rewritten.

Finally, after the two control steps have been executed, the release flatfiles
are generated by running `make_flatfiles.pl`, which creates `MGnifam.seed`,
`MGnifam.hmm` and `MGnifam.hmm.dat` flatfiles.

Moreover, `relnotes`, `README` and `userman.txt` files must be copied from
previous release and `README` file must be updated with the number of families
in the new release.

New release files should be put in `/nfs/production/xfam/mgnifam/releases/`,
under a folder that is named according the new release number.
