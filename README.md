# TrueConsense

TrueConsense is a nucleotide consensus caller capable of creating a biologically valid consensus-sequence in reference-based virus sequencing experiments.

TrueConsense is, in contrast to more common consensus-sequence methodology, not using VCF files. Instead, a consensus-sequence is made directly from a BAM-file in order to get a proper view of positional context.

TrueConsense compensates for common sequencing and/or alignment artefacts, keeps track of open reading frames, and is able to generate consensus-sequences on multiple coverage-thresholds at the same time.

TrueConsense is available under the [AGPLv3 licence](https://www.gnu.org/licenses/agpl-3.0.en.html)