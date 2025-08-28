from TrueConsense.TrueConsense import main


def test_main():
    args = [
        "--input",
        "test/e2e/data/ESIB_EQA_2024_SARS1_01.bam",
        "--output",
        "test/e2e/data/ESIB_EQA_2024_SARS1_01_output.fa",
        "--reference",
        "test/e2e/data/ESIB_EQA_2024_SARS1_01_reference.fasta",
        "--features",
        "test/e2e/data/ESIB_EQA_2024_SARS1_01_features.gff",
        "--coverage-level",
        "30",
        "--samplename",
        "ESIB_EQA_2024_SARS1_01",
    ]

    main(args)

def test2():
    args = [
        "--input",
        "test/e2e/data/ESIB_EQA_2024_SARS1_02.bam",
        "--reference",
        "test/e2e/data/ESIB_EQA_2024_SARS1_02_reference.fasta",
        "--features",
        "test/e2e/data/ESIB_EQA_2024_SARS1_02_features.gff",
        "--coverage-level",
        "30",
        "--samplename",
        "ESIB_EQA_2024_SARS1_02",
        "--output",
        "test/e2e/data/ESIB_EQA_2024_SARS1_02_output.fa",
        "--variants",
        "test/e2e/data/ESIB_EQA_2024_SARS1_02_variants.vcf",
        "--output-gff",
        "test/e2e/data/ESIB_EQA_2024_SARS1_02_output.gff",
        "--depth-of-coverage",
        "test/e2e/data/ESIB_EQA_2024_SARS1_02_depth_of_coverage.tsv",
        "--threads",
        "6"
    ]

    main(args)
