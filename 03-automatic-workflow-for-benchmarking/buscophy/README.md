# Buscophy pipeline for generation of phylogenetic trees from busco genes
All rules and base pipeline is taken from https://gitlab.leibniz-lib.de/smartin/buscophy written by [Sebastian Martin](https://gitlab.leibniz-lib.de/smartin)

## The main changes so far are:

* An update from BUSCO v5 to v6 to be able to incorporate the odb12 databases. 
* Imposing use of --metaeuk to get both nt and aa files
* switch to using mafft for alignment of both aa and nt files

## Things still to be implemented:

* Change the workflow so that individual gene trees are made instead of a tree based on one supermatrix (could be an option)
* Add multi-copy genes and implement a procedure for choosing the paralogs
  - Scripts are now in the `scripts` folder detailing a method for paralog selection based on gene trees
* Add the option to run busco with the `--augustus` argument. There needs to be some change with the augustus config in the singularity image
* Modularise the pipeline to help tidy it up and allow for plug-and-play configuration
* Tidy up the conda envs - at the moment these have no version numbers. They should be consistent with the singularity/docker images.

## In order to run:

Create base conda environment with

`conda env create --file base_snakemake.yml`

I have tested this with conda and singularity. It should also work with apptainer

`snakemake --use-conda --cores 40 #--conda-frontend conda`
`snakemake --cores 40 --use-singularity --singularity-args="--cleanenv --no-home --bind /path/to/cwd/"`

At the moment all parameters should be changed in the `config.yaml` file including the busco lineage, minimum number of species and genes and maximum number of threads per job.

The input genomes should be placed in the `species` folder (can be changed in the `config.yaml`) and should be unzipped and named with the `.fas` extension. For example, as Sebastian suggested for his example set, run the following:

```
mkdir input_test
cd input_test/

# download some lepidoptera genomes
# the original download adress from lepbase.org was not available anymore. This is currently located here:
# https://lepbase.cog.sanger.ac.uk/archive/v4/
wget "https://lepbase.cog.sanger.ac.uk/archive/v4/sequence/Amyelois_transitella_v1_-_scaffolds.fa.gz"
wget "https://lepbase.cog.sanger.ac.uk/archive/v4/sequence/Bicyclus_anynana_v1.2_-_scaffolds.fa.gz"
wget "https://lepbase.cog.sanger.ac.uk/archive/v4/sequence/Bombyx_mori_ASM15162v1_-_scaffolds.fa.gz"
wget "https://lepbase.cog.sanger.ac.uk/archive/v4/sequence/Callimorpha_dominula_k41_-_scaffolds.fa.gz"
wget "https://lepbase.cog.sanger.ac.uk/archive/v4/sequence/Calycopis_cecrops_v1.1_-_scaffolds.fa.gz"
wget "https://lepbase.cog.sanger.ac.uk/archive/v4/sequence/Glyphotaelius_pellucidus_k51_-_scaffolds.fa.gz"

gunzip *fa.gz

#rename genome assemblies
mv Amyelois_transitella_v1_-_scaffolds.fa Amyelois_transitella.fas 
mv Bicyclus_anynana_v1.2_-_scaffolds.fa Bicyclus_anynana.fas
mv Bombyx_mori_ASM15162v1_-_scaffolds.fa Bombyx_mori.fas 
mv Callimorpha_dominula_k41_-_scaffolds.fa Callimorpha_dominula.fas 
mv Calycopis_cecrops_v1.1_-_scaffolds.fa Calycopis_cecrops.fas 
mv Glyphotaelius_pellucidus_k51_-_scaffolds.fa Glyphotaelius_pellucidus.fas 
```
