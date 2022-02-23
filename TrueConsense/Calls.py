import pandas as pd
import pysam
import numpy as np
from .func import chunks, ambig


class Calls:
    ref = None
    p_index = None
    IncludeAmbig = None
    significance = None
    len = None

    def __init__(self, ref, IncludeAmbig=True, significance=0.5):
        self.ref = ref
        self.len = len(ref)
        self.IncludeAmbig = IncludeAmbig
        self.p_index = pd.DataFrame(
            dict(
                refnuc=list(ref),
                query_sequences=np.empty((self.len, 0)).tolist(),
                picked_call=np.nan,
            ),
            index=range(1, self.len + 1),
        )
        self.significance = significance

    def fill_positions_from_bam(self, bamfile, positions=None):
        if positions:
            positions = set(positions)

        bamfile = pysam.AlignmentFile(bamfile, "rb")
        pileup = bamfile.pileup(
            stepper="nofilter", max_depth=10000000, min_base_quality=0
        )
        columns = ["pos", "query_sequences"]
        # 1 Is added to the position because our index starts at 1
        translate_table = str.maketrans("*", "-", "+0123456789")
        updated_p_index = pd.DataFrame(
            (
                (
                    p.pos + 1,
                    [
                        call.split("-")[0].translate(translate_table)
                        for call in p.get_query_sequences(add_indels=True)
                        if call
                    ],
                )
                for p in pileup
                if not positions or p.pos + 1 in positions
            ),
            columns=columns,
        )
        self.p_index = updated_p_index.set_index("pos").combine_first(self.p_index)
        self.p_index["cov"] = self.p_index.query_sequences.apply(len)

    def remove_calls_with_low_coverage(self, mincov):
        self.p_index.at[self.p_index["cov"] < mincov, "calls"] = np.nan

    def insert_calls(self, calls):
        """Inserts calls into p_index['picked_call']

        Args:
            p_index (pd.DataFrame): A pandas dataframe with a column picked_call.
            calls (list(call)): A list of calls that have to be inserted into p_index.

        Returns:
            list(call): The list of calls that were replaced.
        """
        list_of_originals = []
        for call in calls:
            list_of_originals.append(self.p_index.at[call["pos"], "picked_call"])
            self.p_index.at[call["pos"], "picked_call"] = call
        return list_of_originals

    def sort_highest_score(self):
        self.p_index["calls"] = self.p_index.calls.map(
            sort_highest_score, na_action="ignore"
        )

    def calculate_scores(self):
        self.p_index["calls"] = self.p_index[["query_sequences"]].apply(
            lambda l: call_counts(
                l.name,  # l.name is actually the query postition (namely the index of the df)
                l["query_sequences"],
            ),
            axis=1,
        )

    def pick_first_in_calls(self):
        self.p_index["picked_call"] = self.p_index["calls"].map(
            lambda c: c[0] if c and c is not np.nan else None
        )

    def score_coding_sequence(
        self, feature, ends_with_stop_codon=True, starts_with_atg=True
    ):
        """Scores a feature based on the protein it can produce.

        It is: the number of codons in frame / total codons

        Args:
            p_index (pd.DataFrame): A pandas dataframe with a column picked_call.
            feature (pd.Series): A row of a gff as produced with gffpandas dataframe with attributes_to_columns
            ends_with_stop_codon (bool, optional): Whether the last codon of the feature is a stop. Defaults to True.
            starts_with_atg (bool, optional): Whether the first codon of the feature is a start. Defaults to True.

        Returns:
            float: The relative number of codons in-frame.
        """
        # TODO: Add finding of the start codon.
        # TODO: Add finding of alternative stop codon if the original is mutated.
        out_of_frame_counter = 0
        in_frame_counter = 0

        encountered_stop = False
        frame_offset = 0  # If mod. 3 applied, it is the phase of the considered nucleotide relative to the start of the feature. This shifts with insertions and deletions.
        position_in_codon = 0
        last_two_nucs = ""  # The last two nucleotides that have been called. This is needed for (stop)codon recognition

        # TODO: Understand why the -1 is needed here.
        # TODO: Consider sequence before feature.start to find alternative starts.
        # TODO: Consider sequence after feature.end to find alternative stops.
        for i, call in enumerate(
            self.p_index[feature.start - 1 : feature.end]["picked_call"]
        ):

            # Make sure that calls afer regions with uncertainty are considered to be in-frame
            if not call:
                # TODO: fix if there is a framshift (insertion or deletion) after a stretch of N's. An ins or del could shift the frame in step, in stead of out of step.
                frame_offset = 0  # Assume you are back in frame if insuffiecient info
                last_two_nucs = ""
                position_in_codon = (
                    i + 1
                ) % 3  # Reset the position in codon to be in phase
                continue

            seq = call["seq"]
            if seq == "-":
                frame_offset = frame_offset - 1
                continue

            seq_with_prev = last_two_nucs + seq
            last_two_nucs = seq_with_prev[-2:]
            codons = chunks(seq_with_prev[2 - position_in_codon :], 3)
            for codon in codons:
                if len(codon) == 3:
                    if encountered_stop or frame_offset % 3:
                        out_of_frame_counter += 1
                    else:
                        in_frame_counter += 1
                        if codon in ["TAA", "TGA", "TAG"]:
                            encountered_stop = True

            position_in_codon = (position_in_codon + len(seq)) % 3

            if len(seq) > 1:
                frame_offset = frame_offset + len(seq) - 1

        score = in_frame_counter / (in_frame_counter + out_of_frame_counter)
        return score

    def consensus(self):
        cons = self.p_index[["picked_call", "calls"]].apply(
            lambda l: self.call_from_picked(l.picked_call, l.calls),
            axis=1,
        )
        return "".join(cons)

    def call_from_picked(self, picked_call, alt_calls):
        if not picked_call:
            return "N"
        if not self.IncludeAmbig:
            return picked_call["seq"]
        if picked_call["seq"] == "-":
            return "-"
        max_score = max(c["score"] for c in alt_calls)
        ambig_threshold = max_score - 0.1  # Everything within 10% of the most-occuring
        calls = [
            c["seq"]
            for c in alt_calls
            if c["score"] >= ambig_threshold and c["seq"] in ["A", "C", "G", "T"]
        ]
        if len(calls) > 1:
            return ambig(calls)
        return picked_call["seq"]


def sort_highest_score(calls):
    """Sorts calls based on their score property"""
    calls.sort(key=lambda x: x["score"], reverse=True)
    return calls


def call_counts(pos, query_sequences):
    """Produces calls and scores from query sequences.

    The produced calls dictionaries with a position (pos), sequence (seq), read_count (n), score (score), relative score (rel_score) and index (index)

    Args:
        pos (int): The position on the reference
        query_sequences (list(str)): A list of all the sequences found at this position (including insertions and deletions)
        mincov (int): The minimum coverage depth at this position

    Returns:
        list(calls): A list of calls at the position. If cov < mincov returns an empty list.
    """
    cov = len(query_sequences)
    if cov < 1:
        return []

    # This method of counting is faster than the built in collections.Counter class
    queries_without_strand = list(map(str.upper, query_sequences))
    unique_queries = set(queries_without_strand)
    counts = [queries_without_strand.count(q) for q in unique_queries]
    max_count = max(counts)

    # TODO: Add strand bias as a metric.
    # TODO: Consider insertions seperately from nucleotide calls. The first base of the insertion call is not considered in the count of that particular base.

    return [
        dict(
            pos=pos,
            seq=seq,
            n=count,
            score=count / cov,
            rel_score=count / max_count,
            # The relative score is the probability that you would choose the alternative call(s) over the most occuring call(s).
            # e.g. there are 5 calls for the nucleotide C, and 3 for T. The relative score for C is 1.0 and for T 3/5=0.6.
            # The particularities of these relative scores make that they multiply in case of combinations of mutations.
        )
        for seq, count in zip(unique_queries, counts)
    ]
