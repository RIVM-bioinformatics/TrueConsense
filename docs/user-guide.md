# User guide

In order for TrueConsense to keep track of the various aspects of the viral-genome, it is necessary to have enough contextual information to work with.

In order to run TrueConsense you will need the following information/files:

* Your aligned reads in BAM format, preferably created with Minimap2 + samtools
* A reference sequence in FASTA format. This should be the same fasta that you used for your alignment.
* An overview of genomic features matching the **reference** sequence, in GFF3 format
* The samplename


## Basic usage example

In TrueConsense, the `--output`/`-o` flag corresponds to an **output directory** instead of a single file.  
This is because TrueConsense can create multiple consensus-sequences at the same time. These various consensus-sequences will then be generated in the given output directory.  

It's also necessary to provide the samplename to TrueConsense as the output consensus-sequence will follow a format as `{SAMPLENAME}_cov_ge_{COVERAGE_THRESHOLD}.fa` where `{SAMPLENAME}` and `{COVEREAGE_THRESHOLD}` will be replaced by the given information.

When running TrueConsense, it's generally best to provide the same amount of threads as the various coverage-thresholds you wish to generate. For example: if you wish to generate consensus-sequences at thresholds 1, 5, and 10, then it's best to also provide 3 threads.


**Example 1**: Create a consensus-sequence at coverage-thresholds 1, 5, and 10. Using 3 threads.
```bash
trueconsense \
    --input input.bam \
    --reference reference.fasta \
    --features reference-features.gff \
    --samplename example \
    --output test_output/ \
    --coverage-levels 1 5 10 \
    --threads 3
```

This command will generate the following output-files in the `test_output` directory:

```
example_cov_ge_1.fa
example_cov_ge_5.fa
example_cov_ge_10.fa
```

---

Aside from just creating a consensus-sequence, TrueConsense can also provide you with VCF-files, updated GFF files and a coverage overview which you can use for downstream analysis.

??? summary "Generating a coverage overview"
    TrueConsense can generate a per-position depth-of-coverage overview that you can use in a downstream analysis.  
    This file can be generated with the `--depth-of-coverage`/`-doc` flag, the output will be provided in TSV format.

    **Example:**
    ```bash
    trueconsense \
        --input input.bam \
        --reference reference.fasta \
        --features reference-features.gff \
        --samplename example \
        --output test_output/ \
        --coverage-levels 1 5 10 \
        --depth-of-coverage example_coverage.tsv \
        --threads 3
    ```

??? summary "Generating VCF-files per coverage threshold"
    TrueConsense can generate a VCF file per coverage threshold that *matches* the generated consensus-sequence on that particular coverage threshold.  
    With this, you can easily generate the various VCF-files in case a downstream analysis requires VCF input. Or if you wish to share data with other researchers but it's not possible to share (large)Fasta files.

    The VCF-output can be generated with the `--variants`/`-vcf` flag, this flag corresponds to an **output directory** and not a single file. The generated files follow the same naming structure as the normal consensus-sequences.

    **Example:**
    ```bash
    trueconsense \
        --input input.bam \
        --reference reference.fasta \
        --features reference-features.gff \
        --samplename example \
        --output test_output/ \
        --coverage-levels 1 5 10 \
        --variants vcf_output/ \
        --threads 3
    ```

??? summary "Generating updated GFF files per coverage threshold"
    TrueConsense keeps track of the open reading frames given in the input-gff and determines new stop-positions and/or new start-positions when applicable for every open reading frame in the input.  
    This is done to make sure the generated consensus-sequence is possible when it comes to virus-biology.  

    The updated GFF-files can be generated with the `--output-gff`/`-ogff` flag, this flag corresponds to an **output directory** and not a single file. The generated files follow the same naming structure as the normal consensus-sequences.

    **Example:**
    ```bash
    trueconsense \
        --input input.bam \
        --reference reference.fasta \
        --features reference-features.gff \
        --samplename example \
        --output test_output/ \
        --coverage-levels 1 5 10 \
        --output-gff gff_output/ \
        --threads 3
    ```

TrueConsense calls IUPAC nucleotide ambiguity-codes by default when an aligned-position has an even (or near-even) split of nucleotides.  
This can be turned off by providing the `--noambiguity`/`-noambig` flag. Please note that this will cause TrueConsense to choose the most prominent nucleotide on a split-position, if the split is *exactly* even on such a position then a random choice will be made between the two (or three) possibilities.

---

## Limitations

TrueConsense only works with alignments in BAM-format as an input. Other inputs such as SAM or CRAM are currently not supported.

Additionally, the generated VCF-files are formatted in a way that they can be used to accurately reconstruct the consensus-sequence when sharing the VCF files with other researchers/institutes/etc with a tool such as `bcftools consensus`. That also means that these VCF-files may not *always* follow the specifications for VCF-files, for example when ambiguity nucleotides are present in the VCF-file.  
The generated VCF-files are made in a way that a tool such as `bcftools` can use it, but these VCF-files may not always work in other tools such as IGV.

TrueConsense can only compensate for common alignment artefacts to a certain degree and does so based on the index that TrueConsense makes which contains a distribution of  "nucleotide events" per column of the alignment.  
There are some edge-cases where the alignment cannot be done without introducing an artefact which cannot be compensated by TrueConsense as information will not be included by the aligner on certain positions.  
In this case, the index on these positions can be overwritten with an index made in an earlier stage of analysis.  
This can be done with the `--index-override` flag, the data must be provided as a compressed CSV (.csv.gz) file. Below is an example of the necessary format, where "X" is deletions and "I" is number of reads with an insertion. The first column (position) has no header.  

||coverage                     |A     |T                                            |C  |G  |X  |I  |
|------|-----------------------|------|---------------------------------------------|---|---|---|---|
|1     |1                      |2     |7                                            |3  |5  |1  |3  |
|2     |1                      |2     |3                                            |4  |10 |0  |0  |

!!! warning "Please only use this when absolutely necessary"
    Using the index override may solve a very specific issue for you particular analysis, but it will also cause the result to be much harder to validate.  
    Additionally, this may introduce new issues as there is little to no validation of the override-data which will be handled as "truthful" data.