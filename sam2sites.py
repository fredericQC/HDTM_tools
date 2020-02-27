#!/usr/bin/env python
# coding: utf-8

from path import Path
import numpy as np
import pandas as pd
import click


def parse_sam_file(samfile, threshold=30) -> pd.DataFrame:
    # import the data
    df = pd.read_csv(
        samfile,
        dtype={
            "chr_name": str,
            "pos": "Int64",
            "flag": "Int64",
            "name": str,
            "score": "Int64",
            "cigar": str
        },
        sep="\t",
        header=None,
        usecols=range(6),
        error_bad_lines=False,
        names=["name", "flag", "chr_name", "pos", "score", "cigar"],
    )
    # get only stranded values
    df = df[df.flag.isin([0, 16])]

    # map flag 0 to strand "+", flag 16 to strand "-"
    df["strand"] = df.flag.replace(
        {0: pd.Categorical("+", ["+", "-"]), 16: pd.Categorical("-", ["+", "-"])}
    )
    # use regex to extract each int/marker and build a secondary df holding these info.
    #
    # typic input:
    # ------------
    #
    # 20      42M
    # 25    72M3S
    # Name: cigar
    #
    # typic output:
    # -------------
    #
    #           size marker
    #    match
    # 20 0        42      M
    # 25 0        72      M
    #    1         3      S
    cigar_df = df.cigar.str.extractall(r"(?P<size>\d+)(?P<marker>[A-Z]{1})")
    cigar_df["size"] = cigar_df["size"].astype(int)

    # because needed to test validity, extract first and last cigar marker.
    df["cigar_first_marker"] = cigar_df.marker.groupby(level=0).first()[df.index]
    df["cigar_last_marker"] = cigar_df.marker.groupby(level=0).last()[df.index]

    # the size can be computed with the sum of the size associated with the
    # markers 'M', 'D', 'N', '=' and 'X'

    df["aligned_read_length"] = (
        cigar_df.loc[cigar_df.marker.str.contains("[MDN=X]")].groupby(level=0).sum()
    )
    # two case of validity: positive strand means we need first 'M' (match) marker,
    # negative strand means we need last 'M' (match) marker
    is_valid_positive = (df.strand == "+") & (df.cigar_first_marker == "M")
    is_valid_negative = (df.strand == "-") & (df.cigar_last_marker == "M")
    # we accept reads that are either a valid positive, negative positive and which
    # length are more than a threshold.
    df["valid"] = (is_valid_positive | is_valid_negative) & (
        df.aligned_read_length > threshold
    )

    valid_df = df[df.valid]
    valid_df["start"] = np.where(
        valid_df.strand == "+",
        valid_df.pos + 3,
        valid_df.pos + valid_df.aligned_read_length - 6,
    )
    return valid_df


def sam_to_site(sam_df: pd.DataFrame) -> pd.DataFrame:
    # we aggregate the scores by start
    score_pos = (
        sam_df[sam_df.strand == "+"].groupby("start").score.count().rename("score_pos")
    )
    score_neg = (
        sam_df[sam_df.strand == "-"].groupby("start").score.count().rename("score_neg")
    )
    score_total = sam_df.groupby("start").score.count().rename("score_total")
    # we take the first chr_name. Should be the same for each occurence
    chr_name = sam_df.groupby("start").chr_name.first()

    # Concatenate the aggregates
    site_df = pd.concat(
        [chr_name, score_pos, score_neg, score_total], axis=1
    ).reset_index()
    # fill the last values
    site_df["end"] = site_df.start + 1
    site_df["name"] = "i"
    # set the empry score to 0
    site_df = site_df.fillna({"score_pos": 0, "score_neg": 0, "score_total": 0})
    # reorder the columns and sort by start values
    site_df = site_df[
        ["chr_name", "start", "end", "name", "score_pos", "score_neg"]
    ].sort_values(by=["chr_name","start"])
    return site_df


def site_to_interval(
    site_df: pd.DataFrame, *, normalization_value: int, score_threshold=1
) -> pd.DataFrame:
    total_reads = site_df["score_pos"].sum() + site_df["score_neg"].sum()

    pos_interval_df = (
        site_df[site_df["score_pos"] > score_threshold]
        .drop("score_neg", axis=1)
        .rename(columns={"score_pos": "score"})
    )
    pos_interval_df["strand"] = "+"

    neg_interval_df = (
        site_df[site_df["score_neg"] > score_threshold]
        .drop("score_pos", axis=1)
        .rename(columns={"score_neg": "score"})
    )
    neg_interval_df["strand"] = "-"

    interval_df = pd.concat([pos_interval_df, neg_interval_df])
    interval_df["normalized_score"] = (
        (interval_df.score / interval_df.score.sum() * normalization_value)
        .astype(int)
        .clip(1)
    )
    interval_df["fake_score"] = 999
    interval_df = interval_df.sort_values(by=["chr_name","start"])
    return interval_df


def build_bed_file(interval_df: pd.DataFrame, filename: Path):
    """
    bed_file
    --------

    "chr_name", "start", "end", "name", "fake_score", "strand"

    pVCR94deltaX3deltaacr2FRT	26	27	i	999	+
    pVCR94deltaX3deltaacr2FRT	47	48	i	999	+
    pVCR94deltaX3deltaacr2FRT	68	69	i	999	+
    pVCR94deltaX3deltaacr2FRT	97	98	i	999	-
    pVCR94deltaX3deltaacr2FRT	245	246	i	999	-
    """
    interval_df[["chr_name", "start", "end", "name", "fake_score", "strand"]].to_csv(
        filename, header=False, index=False, sep="\t"
    )


def build_stranded_unnormalized_bedgraph_file(
    interval_df: pd.DataFrame, filename: Path
):
    """
    stranded_unnormalized_bedgraph_file
    --------

    "chr_name", "start", "end", "score"

    pVCR94deltaX3deltaacr2FRT	26	27	2400
    pVCR94deltaX3deltaacr2FRT	47	48	60
    pVCR94deltaX3deltaacr2FRT	68	69	60
    pVCR94deltaX3deltaacr2FRT	97	98	1380
    pVCR94deltaX3deltaacr2FRT	245	246	3780
    """
    interval_df[["chr_name", "start", "end", "score"]].to_csv(
        filename, header=False, index=False, sep="\t"
    )


def build_stranded_bedgraph_file(interval_df: pd.DataFrame, filename: Path):
    """
    stranded_bedgraph_file
    --------

    "chr_name", "start", "end", "normalized_score"

    pVCR94deltaX3deltaacr2FRT	26	27	2189
    pVCR94deltaX3deltaacr2FRT	47	48	54
    pVCR94deltaX3deltaacr2FRT	47	48	54
    pVCR94deltaX3deltaacr2FRT	68	69	54
    pVCR94deltaX3deltaacr2FRT	97	98	1258
    pVCR94deltaX3deltaacr2FRT	245	246	3448
    """
    interval_df[["chr_name", "start", "end", "normalized_score"]].to_csv(
        filename, header=False, index=False, sep="\t"
    )

def build_unstranded_bedgraph_file(interval_df: pd.DataFrame, filename: Path):
    """
    unstranded_bedgraph_file
    --------

    "chr_name", "start", "end", "normalized_score"

    pVCR94deltaX3deltaacr2FRT	26	27	2189
    pVCR94deltaX3deltaacr2FRT	47	48	108<-(score for insertions on the same position, but different strands are merged)
    pVCR94deltaX3deltaacr2FRT	68	69	54
    pVCR94deltaX3deltaacr2FRT	97	98	1258
    pVCR94deltaX3deltaacr2FRT	245	246	3448
    """

    interval_df = interval_df.groupby(["chr_name", "start", "end"],sort=False).agg({'normalized_score':'sum'}).reset_index()
    interval_df[["chr_name", "start", "end", "normalized_score"]].to_csv(
        filename, header=False, index=False, sep="\t"
    )


def build_scored_bed_file(interval_df: pd.DataFrame, filename: Path):
    """
    scored_bed_file
    --------

    "chr_name", "start", "end", "name", "score", "strand"

    pVCR94deltaX3deltaacr2FRT	26	27	i	2400	+
    pVCR94deltaX3deltaacr2FRT	47	48	i	60	+
    pVCR94deltaX3deltaacr2FRT	68	69	i	60	+
    pVCR94deltaX3deltaacr2FRT	97	98	i	1380	-
    pVCR94deltaX3deltaacr2FRT	245	246	i	3780	-
    """
    interval_df[["chr_name", "start", "end", "name", "normalized_score", "strand"]].to_csv(
        filename, header=False, index=False, sep="\t"
    )


@click.command()
@click.argument("input_filename", type=click.File("r"))
@click.argument("normalization_value", type=int)
@click.option("--output_dir", type=click.Path(exists=True), default=Path("."))
@click.option("--read_len_threshold", default=30, type=int)
@click.option("--score_threshold", default=1, type=float)
@click.option("--create_bedfile", default=True)
@click.option("--create_scored_bedfile", default=True)
@click.option("--create_stranded_bedgraph", default=True)
@click.option("--create_unstranded_bedgraph", default=True)
@click.option("--create_stranded_unnormalized_bedgraph", default=True)
def cli_interface(
    input_filename: Path,
    normalization_value: int,
    output_dir: Path = Path("."),
    read_len_threshold: int = 30,
    score_threshold: float = 1,
    create_bedfile=True,
    create_scored_bedfile=True,
    create_stranded_bedgraph=True,
    create_unstranded_bedgraph=True,
    create_stranded_unnormalized_bedgraph=True,
):
    output_dir = Path(output_dir)
    basename = Path(input_filename.name).basename().stripext()
    if basename == "<stdin>":
        basename = "cigar_out"
    sam_df = parse_sam_file(input_filename, read_len_threshold)
    site_df = sam_to_site(sam_df)
    interval_df = site_to_interval(
        site_df,
        normalization_value=normalization_value,
        score_threshold=score_threshold,
    )
    if create_bedfile:
        build_bed_file(interval_df, output_dir / "%s.bed" % basename)
    if create_scored_bedfile:
        build_scored_bed_file(interval_df, output_dir / "%s_scored.bed" % basename)
    if create_stranded_bedgraph:
        build_stranded_bedgraph_file(
            interval_df, output_dir / "%s_stranded.bg" % basename
        )
    if create_unstranded_bedgraph:
        build_unstranded_bedgraph_file(
            interval_df, output_dir / "%s_unstranded.bg" % basename
        )
    if create_stranded_unnormalized_bedgraph:
        build_stranded_unnormalized_bedgraph_file(
            interval_df, output_dir / "%s_stranded_unnormalized.bg" % basename
        )


if __name__ == "__main__":
    cli_interface()
