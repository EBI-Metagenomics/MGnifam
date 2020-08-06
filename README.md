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

## Setup and usage

### MGnifam setup (only for old release)

First step for setting up MGnifam database is by downloading the Pfam source
code directly from GitHub as follows:

```
git clone https://github.com/ProteinsWebTeam/Pfam.git
cd Pfam
git fetch origin mgnifam
git checkout mgnifam
```

MGnifam core files can be found at:
  1. `Pfam/PfamScripts/mgnifam/` for building MGnifam clusters;
  2. `Pfam/PfamScripts/mgnifam_release` for releasing MGnifam clusters;
  3. `nfs/production/xfam/mgnifam/Conf/mgnifam.conf` (example path) for MGnifam configuration;

### Environmental setup

Both the old and the new MGnifam release pipelines must call external scripts:
while the former is very dependent from the MGnifam Perl scripts (such as
mgseed, pfbuild, mfbuild, ...) and its configuration (mgnifam.conf), the latter
depends on only a few external scripts such as MobiDB Lite and HMMER3.

In both cases, however, it could be that scripts require specific environmental
variables setup to run properly (e.g. MobiDB Lite requires both Python 2 and 3
to avoid crashing). For that purpose, boh pipelines implement the
`[-e ENV_PATH]` optional parameter, allowing users to define a JSON file
containing a dictionary: environmental variables names are the keys, while
values are lists of strings which will the be automatically concatenated by the
command line script using `: (double dots)` as separator.

Example of `env.json` file:

```
{
  "PATH": ["first/path", "second/path", "third/path"],
  "PYTHONPATH": ["first/python/path", "second/python/path"]
}
```

### Input

Input of the pipelines are one or multiple file holding table in it, whose first
column is cluster names (and second is cluster size for example, if any). Then,
a batch size can be specified (how many HMM build must be run in parallel) and
a maximum number of clusters can be set, limiting the number of clusters to
read from input files (sequentially).

### Annotator

Each pipeline requires an author name: it is required to identify who
actually ran the pipeline (either the new or the old one) and set it into
clusters DESC files.

### Version 1.0.0 (The old one)

Old version of MGnifam release pipeline is strongly dependent to MGnifam source
Perl scripts. Because of that, reference to them has to be kept inside of
environmental variables, as well as a reference to the MGnifam configuration
file, which is used to retrieve all the settings (such as database credentials,
paths to database, ...).

**Usage**

```
release100.py [-h] --in_path IN_PATH [IN_PATH ...] --out_path OUT_PATH
                   -a AUTHOR_NAME [-n MAX_CLUSTERS] [-b BATCH_SIZE]
                   [-v VERBOSE] [-e ENV_PATH]
```

### Version 2.0.0 (The new one)

New version of MGnifam release pipeline is not dependent on MGnifam Perl
script. However, it still uses external scripts under the hood, while much less
intensively than the old pipeline, such as MobiDB Lite and HMMER3. The former
one, specifically, requires careful setup of both Python 2 and Python 3
dependencies in order to run multiple disorder predictors being it just an
ensemble predictor, whose inner predictor have different requirements.

In this pipeline, batch size has a very high importance: in fact, parallel
computations results are all sent back to main process, filling its memory.
For example, one of the most memory intensive step is fasta sequences search
against UniProt or MGnifam datasets, which will retrieve entire fasta sequences
in dictionaries. In order to avoid memory issues is therefore mandatory to
allocate a reasonable amount of memory to the main process and limiting the
number of parallel computations (batch size) to a number between 10 and 100
thousands.

**Requirements**

- python >= 3.7.x
- dask >= 2.17.2
- dask_jobqueue >= 0.7.1
- numpy >= 1.18.1
- muscle >= 3.8.31
- mobidb_lite >= 2.0.x
- hmmer >= 3.1b2
  - hmmbuild
  - hmmsearch
  - hmmalign

Note: MobiDB Lite script reequires both Python 2.x and 3.x. While pipeline
can be run directly from Python3.x, user must provide a suitable environment to
subprocesses (passing a properly formatted JSON file to
`--environ_path [ENVIRON_PATH]` option).

**Usage**

```
release200.py [-h] --in_path IN_PATH [IN_PATH ...] --out_path OUT_PATH
                   -a AUTHOR_NAME [-n MAX_CLUSTERS] [-b BATCH_SIZE]
                   [-v VERBOSE] [-e ENV_PATH]
```
