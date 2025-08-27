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
