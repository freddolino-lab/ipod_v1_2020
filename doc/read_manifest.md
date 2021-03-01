A read manifest file (`read_manifest.txt`) must be present in each of the IPOD-HR data directories, and contains information on all of the raw read files that should be processed for that condition/sample-type combination. The file is a whitespace-delimited plain text file that must contain five fields:

1. Forward read file name (path relative to the data directory). This can be a fastq file, optionally compressed via bzip2 or gzip; if compressed, the file name should end in `.gz` or `.bz2`, as appropriate.
2. Reverse read file name (reads should be in the same order as the forward read file). **Note that IPOD-HR requires paired-end reads**
3. Adapter sequence to be trimmed from the 3' end of reads. The truseq shared adapter sequence (`AGATCGGAAGAGC`) often suffices.
4. PHRED score base. Should almost always be 33 or 64
5. Sample name. This name will be applied to all output files produced from this sample

It is good, although not required, practice to put all raw read files in a `raw/` subdirectory of the data directory, in which case the forward and reverse read files listed in the manifest should begin with `raw/` to reflect their path relative to the data directory.

