The sample-level configuration files are a key component of the overall layout of the IPOD-HR pipeline. These files are structured as [TOML](https://toml.io/en/) documents, and contain the sections documented below. The sample dataset distributed with the IPOD-HR pipeline provides a practical example of what these files should contain.

# general
This section contains overall information on the sample and input/output locations.

`bs_dirname`: Directory name for read processing and bootstrapping
`bs_suffix`: Suffix for bootstrap-resampled reads
`orig_suffix`: Suffix for original (non-bootstrapped) data
`out_prefix`: Sample-specific prefix to be applied to all processed output files for this sample
`output_path`: Subdirectory of the sample directory to write processed data to
`pretty_name`: Human-readable description of this sample

# bootstrap
This section contains information on the bootstrap resampling that can be used for some error estimation.

`samtools_threads`: Number of threads to use while processing the aligned reads with `samtools`
`bootstrap_threads`: Number of threads to use during bootstrap resampling of reads

# ipod
This section contains information on the data files in the data directory for the IPOD-HR samples. 

`directory`: The name of the directory containing the IPOD-HR samples; this must be a subdirectory of the current sample directory.
`sample_names`: Names of the IPOD-HR samples to be included in the analysis; these must correspond to **sample names** in the read manifest for the IPOD-HR files.

# inp
This section is structed as for the `ipod` section above, but for the input samples.

# chip
This section is structed as for the `ipod` section above, but for the RNA polymerase ChIP-seq samples.
