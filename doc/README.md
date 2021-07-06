% IPOD-HR Analysis Pipeline 

# Introduction

This source tree contains a complete system suitable for running read processing and scoring on IPOD data sets, using version 1.2 of the processing pipeline (currently in revision; an outgrowth of the methods described in [this paper](https://doi.org/10.1101/2020.01.29.924811)). Note that several hardcoded choices in the present version make it suitable only for analysis of data from the *E. coli* MG1655 strain, using the U00096.3 reference genome. Ongoing development is in progress for a more general purpose version of these workflows.

# Installation

The analysis pipeline provided here is reliant on several excellent pieces of open source software, and in some cases requires specific versions in order to function properly. To simplify the process of establishing a compatible environment, we provide a [conda](https://docs.conda.io/en/latest/) environment definition in  the accompanying file `ipod_conda_environment.yml` that will provide nearly all tools necessary for running IPOD-HR; this can be used as a checklist (and guide to specific required versions) even if not working in a conda environment. In addition to the software noted there, [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) must be installed and in the active PATH; all development, testing, and application of the pipeline was performed using bowtie2 version 2.1.0.

As a simplified alternative for getting a functional environment for analyzing IPOD-HR data, we also provide a [singularity](https://sylabs.io/guides/3.7/user-guide/) container that has a fully functional IPOD-HR analysis environment. Instructions for using the containerized environment are given [below](#containerized-version).

# File preparation

The IPOD processing pipeline requires a strict file hierarchy in order to properly find and process the data files.

## Top level directory
At the top level directory, there must be a text file containing directory/configuration pairs, formatted as
    directory1 config1
    directory2 config2
    directory3 config3
etc...

All of the specified directories must exist in the top-level directory, and must contain a configuration file corresponding to the config name given. We refer to each of the specified subdirectories as *sample directories* in the documentation below; typically each one will represent all of the data present on a single biological condition.

## Sample directories

Each sample directory must contain a sample-level configuration file (as noted above), as well as three subdirectories named `ipod`, `inp`, and `chip`. We collectively refer to those subdirectories as the *data directories*; they correspond to the IPOD-HR, input, and RNA polymerase ChIP-seq samples, respectively, under the biological condition being considered.

### Sample-level configuration file

The sample-level configuration file is a [TOML](https://toml.io/en/) document; formatting details are given in the accompanying sample_config.md file.

### Data directories

Each data directory must contain:
* A read manifest file named `read_manifest.txt` (documented in read_manifest.md)
* Directories labeled `raw`, `aligned`, and `bootstrap`
* Within the `raw` directory, the read files reference from the read manifest

# Running the analysis

Assuming that the installation prerequisites described above have been met, and that this source code distribution is present in a directory called `{SRC_LOC}`, the entire pipeline can be run using the python program in `{SRC_LOC}/drivers/run_all.py`, specifying as a single command line argument the name of the top-level directory listing noted above. 

# Output files

The IPOD-HR analysis pipeline will produce several intermediate files as well as a final set of results. Intermediate files are typically the results of individual pipeline steps (e.g., running `bowtie`). The final results will be written to the directory specified in `general -> out_prefix` option of the sample-level configuration file. Within that directory, the files of typical interest are:

* `OUTPREFIX_v6rz_chipsub.gr` -- .gr file containing the robust Z-scores after ChIP subtraction. This is the most commonly used output in practice.
* `OUTPREFIX_v6rzlog10p_chipsub.gr` -- .gr file containing the log10p-scaled robust Z-scores, effectively yielding p-values assuming a standard normal null distribution
* `OUTPREFIX_ipod_vs_inp_lograt.gr` -- .gr file of the IPOD/input log ratios, prior to subtraction of RNA polymerase occupancy
* `OUTPREFIX_chip_vs_inp_lograt.gr` -- .gr file of the RNA Polymerase ChIP/input log ratios


# Example data

To provide a test case that illustrates the full set of inputs needed to apply the IPOD-HR analysis tools in a nontrivial case, we distribute a bundle referred to as `IPODHR_MINIMAL_TEST.tbz` which contains the complete directory structure and raw reads needed to process an IPOD-HR experiment; in this case, all data are taken from WT MG1655 growing in log phase in either rich or minimal media. The full test data set can be downloaded [here](https://drive.google.com/file/d/13G8r3cPTBloMF3Bl-U5T_jwrHgsb4GxN/view?usp=sharing). Users will find examples for all required configuration/input files, and can also run the complete analysis of the test data set by entering the sample data directory and calling

`python {SRC_PATH}/drivers/run_all.py all_conditions.txt`

Where {SRC_PATH} indicates the location of the analysis code distributed here.

Included in our sample data distribution are gold standard files for the final results generated by running the pipeline, obtained using our development environment. The results obtained from a test run by calling
`make diff`
which will display the magnitude of the differences between the files generated in your run, and those present in the gold standard files. Both the RMSD and MAX differences **should** be zero if the software environment has been appropriately reproduced.

Note that while the data sets used in this test case are relatively small, they still are designed to provide a non-trivial working example, and will likely take several hours to run on a decently powerful workstation.

# Containerized version

As an alternative to allow rapid and reproducible setup of the IPOD-HR postprocessing pipeline described here, we have also made a [singularity container](https://drive.google.com/file/d/1ruiB1IjLUQxaZjNn4GXVRZmQ08wDntsO/view?usp=sharing) available that provides a complete, self-enclosed environment for data analysis. The environment can be entered by calling `singularity run ipod_v1.2.sif`; the IPOD-HR analysis source code tree will then be mounted at `/src_for_distrib_dec2020`. We highly recommend familiarizing yourself with fundamental concepts in singularity containers prior to using this environment; in particular, it is likely necessary to mount the directory containing your data tree so that it is accessible within the container. As an example session, running on an Ubuntu 14.04 host operating system with singularity version 3.7.1 as the guest, the singularity container could be invoked with

`singularity run -B /data/petefred/TEST_MINIMAL:/testdata -B /run/shm:/run/shm ipod_v1.2.sif`

Here `/data/petefred/TEST_MINIMAL` is the location of the test data described above, which will then be mounted at `/testdata` in the container environment. The `-B /run/shm:/run/shm` motif is necessary to permit proper functioning of the python `multiprocessing` module on an Ubuntu host, and may not be necessary in other environments (e.g., we have not found it to be needed on a Red Hat host OS). 

Having entered the singularity environment, the test case could then be run by executing

`cd /testdata; python /src_for_distrib_dec2020/drivers/run_all.py all_conditions.txt > test.log 2> test.err`

After the run is complete, the results should be checked using the procedure described [above](#example-data). 

# Postprocessing tools

We also include in this source code distribution the python programs needed for key postprocessing tasks from the accompamying manuscripts, namely those used for calling individual TF binding peaks, calling extended protein occupancy domains (EPODs), and for plotting and consensus clustering of TFs based on their binding profiles. Documentation for these programs is included in the postprocessing.md file in this directory. 
