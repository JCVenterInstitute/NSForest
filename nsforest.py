#!/usr/bin/env python

import argparse
import os
from pathlib import Path

import scanpy as sc

import nsforest as ns
from nsforest import nsforesting


def run_nsforest_with_preprocessing(
    h5ad_filepath,
    cluster_header="cell_type",
    cluster_list=[],
    results_dirpath=".",
    total_counts=0,
):
    """Run NSForest using the specified dataset filepath, and
    cluster_header first preprocessing the dataset.

    Parameters
    ----------
    h5ad_filepath : str
        The dataset filepath
    cluster_header : str
        The cluster header
    cluster_list: list
        List of clusters, default: all clusters
    results_dirpath : str
        Directory path for the results
    total_counts : int
        Total counts after downsampling, default: 0, do not downsample

    Returns
    -------
    None

    Notes
    -----
    - Some datasets have multiple annotations per sample
    (ex. "broad_cell_type" and "granular_cell_type"). NSForest can be
    run on multiple `cluster_header`'s. Combining the parent and child
    markers may improve classification results.

    - `adata.var_names` must be unique. If there is a problem, usually
    it can be solved by assigning `adata.var.index =
    adata.var["ensembl_id"]`.

    - Some datasets are too large and need to be downsampled to be run
    through the pipeline. When downsampling, be sure to have all the
    granular cluster annotations represented.

    - Only run ns.pp.dendrogram() if there is no pre-defined
    dendrogram order. This step can still be run with no effects, but
    the runtime may increase.
    """
    # Assign results filename and directory
    h5ad_filename = os.path.basename(h5ad_filepath)
    pp_h5ad_filename = prefix_extension(h5ad_filename, "_pp")
    pp_h5ad_filepath = f"{results_dirpath}/{pp_h5ad_filename}"

    # Run NSForest if results do not exist
    if not os.path.exists(pp_h5ad_filepath):

        print(f"Loading unprocessed AnnData file: {h5ad_filename}")
        up_adata = sc.read_h5ad(h5ad_filepath)

        # Do not downsample by default
        ds_adata = up_adata  # No need to copy
        if total_counts > 0:

            print("Calculating QC metrics")
            up_metrics = sc.pp.calculate_qc_metrics(up_adata)
            if up_metrics[1]["total_counts"].sum() > total_counts:

                print("Downsampling unprocessed AnnData file")
                ds_adata = sc.pp.downsample_counts(
                    up_adata, total_counts=total_counts, copy=True
                )

        print("Generating scanpy dendrogram")
        # Dendrogram order is stored in
        # `pp_adata.uns["dendrogram_cluster"]["categories_ordered"]`
        pp_adata = ds_adata.copy()
        pp_adata.obs[cluster_header] = pp_adata.obs[cluster_header].astype(str)
        pp_adata.obs[cluster_header] = pp_adata.obs[cluster_header].astype("category")
        pp_adata = ns.pp.dendrogram(
            pp_adata,
            cluster_header,
            save=False,
            output_folder=results_dirpath,
            outputfilename_suffix=cluster_header,
        )

        print("Calculating cluster medians per gene")
        pp_adata = ns.pp.prep_medians(pp_adata, cluster_header)

        print("Calculating binary scores per gene per cluster")
        pp_adata = ns.pp.prep_binary_scores(pp_adata, cluster_header)

        print(f"Saving preprocessed AnnData file: {pp_h5ad_filename}")
        pp_adata.write_h5ad(pp_h5ad_filepath)

        print(f"Running NSForest for preprocessed AnnData file: {pp_h5ad_filename}")

        results = nsforesting.NSForest(
            pp_adata,
            cluster_header,
            cluster_list=cluster_list if isinstance(cluster_list, list) else [cluster_list],
            output_folder=f"{results_dirpath}/",
            outputfilename_prefix=cluster_header,
        )

    else:
        print(f"Completed NSForest for preprocessed AnnData file: {pp_h5ad_filename}")


def preprocess_adata_file(
    inp_adata_file, cluster_header, out_adata_file, results_dirpath="."
):
    """Preprocess the specified dataset filepath.

    Parameters
    ----------
    inp_adata_file : str
        Input AnnData file name
    cluster_header : str
        The cluster header
    out_adata_file : str
        Output AnnData file name
    results_dirpath : str
        Directory path for the results

    Returns
    -------
    None

    """
    print(f"Loading input AnnData file: {inp_adata_file}")
    inp_adata = sc.read_h5ad(inp_adata_file)

    print("Generating scanpy dendrogram")
    # Dendrogram order is stored in
    # `out_adata.uns["dendrogram_cluster"]["categories_ordered"]`
    out_adata = inp_adata.copy()
    out_adata.obs[cluster_header] = out_adata.obs[cluster_header].astype(str)
    out_adata.obs[cluster_header] = out_adata.obs[cluster_header].astype("category")
    out_adata = ns.pp.dendrogram(
        out_adata,
        cluster_header,
        save=False,
        output_folder=results_dirpath,
        outputfilename_suffix=cluster_header,
    )

    print("Calculating cluster medians per gene")
    out_adata = ns.pp.prep_medians(out_adata, cluster_header)

    print("Calculating binary scores per gene per cluster")
    out_adata = ns.pp.prep_binary_scores(out_adata, cluster_header)

    print("Writing clusters")
    with open("clusters.txt", "w") as fp:
        for cluster in out_adata.obs[cluster_header].cat.categories:
            fp.write(f"{cluster}\n")

    print(f"Saving output AnnData file: {out_adata_file}")
    out_adata.write_h5ad(out_adata_file)


def downsample_adata_file(inp_adata_file, total_counts, out_adata_file):
    """Calculates quality control metrics, and sums total counts for
    each cell to downsample, if needed.

    Parameters
    ----------
    inp_adata_file : str
        Input AnnData file name
    total_counts : int
        Maximum sum of total counts for each cell
    out_adata_file : str
        Output AnnData file name

    Returns
    -------
    None

    See:
        https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html
    """
    print(f"Loading input AnnData file: {inp_adata_file}")
    inp_adata = sc.read_h5ad(inp_adata_file)

    print("Calculating QC metrics")
    metrics = sc.pp.calculate_qc_metrics(inp_adata)
    if metrics[1]["total_counts"].sum() > total_counts:
        print("Downsampling AnnData file")
        out_adata = sc.pp.downsample_counts(
            inp_adata, total_counts=total_counts, copy=True
        )

    else:
        out_adata = inp_adata  # No need to copy

    print(f"Saving output AnnData file: {out_adata_file}")
    out_adata.write_h5ad(out_adata_file)


def generate_scanpy_dendrogram(
    inp_adata_file, cluster_header, out_adata_file, results_dirpath="."
):
    """Plots a dendrogram of the categories defined in grouping by
    cluster_header.

    Parameters
    ----------
    inp_adata_file : str
        Input AnnData file name
    cluster_header : str
        Column in `adata.obs` storing cell annotation
    results_dirpath : str
        Directory path for the results
    out_adata_file : str
        Output AnnData file name

    Returns
    -------
    None

    See:
        https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.dendrogram.html
    """
    print(f"Loading input AnnData file: {inp_adata_file}")
    inp_adata = sc.read_h5ad(inp_adata_file)

    print("Generating scanpy dendrogram")
    # Dendrogram order is stored in:
    #   `out_adata.uns["dendrogram_cluster"]["categories_ordered"]`
    out_adata = inp_adata.copy()
    out_adata.obs[cluster_header] = out_adata.obs[cluster_header].astype(str)
    out_adata.obs[cluster_header] = out_adata.obs[cluster_header].astype("category")
    out_adata = ns.pp.dendrogram(
        out_adata,
        cluster_header,
        save=False,
        output_folder=results_dirpath,
        outputfilename_suffix=cluster_header,
    )

    print(f"Saving output AnnData file: {out_adata_file}")
    out_adata.write_h5ad(out_adata_file)


def calculate_cluster_medians_per_gene(inp_adata_file, cluster_header, out_adata_file):
    """Calculate the median expression matrix.

    Parameters
    ----------
    inp_adata_file : str
        Input AnnData file name
    cluster_header : str
        Column in `adata.obs` storing cell annotation
    out_adata_file : str
        Output AnnData file name

    Returns
    -------
    None

    See:
        https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.dendrogram.html
    """
    print(f"Loading input AnnData file: {inp_adata_file}")
    inp_adata = sc.read_h5ad(inp_adata_file)

    print("Calculating cluster medians per gene")
    out_adata = ns.pp.prep_medians(inp_adata, cluster_header)

    print(f"Saving output AnnData file: {out_adata_file}")
    out_adata.write_h5ad(out_adata_file)


def calculate_binary_scores_per_gene_per_cluster(
    inp_adata_file, cluster_header, out_adata_file
):
    """Calculating the binary scores of each gene per
    `cluster_header`.

    Parameters
    ----------
    inp_adata_file : str
        Input AnnData file name
    cluster_header : str
        Column in `adata.obs` storing cell annotation
    out_adata_file : str
        Output AnnData file name

    Returns
    -------
    None
    """
    print(f"Loading input AnnData file: {inp_adata_file}")
    inp_adata = sc.read_h5ad(inp_adata_file)

    print("Calculating binary scores per gene per cluster")
    out_adata = ns.pp.prep_binary_scores(inp_adata, cluster_header)

    print(f"Saving output AnnData file: {out_adata_file}")
    out_adata.write_h5ad(out_adata_file)


def run_nsforest_without_preprocessing(
        inp_adata_file, cluster_header, cluster_list=[], results_dirpath="."
):
    """Performs the main NS-Forest algorithm to find a list of
    NS-Forest markers for each `cluster_header`.

    Parameters
    ----------
    inp_adata_file : str
        Input AnnData file name
    cluster_header : str
        Column in `adata.obs` storing cell annotation
    results_dirpath : str
        Directory path for the results
    cluster_list: list
        List of clusters, default: all clusters

    Returns
    -------
    None
    """
    print(f"Loading input AnnData file: {inp_adata_file}")
    inp_adata = sc.read_h5ad(inp_adata_file)

    print(f"Running NSForest for preprocessed AnnData file: {inp_adata_file}")
    nsforesting.NSForest(
        inp_adata,
        cluster_header,
        cluster_list=cluster_list if isinstance(cluster_list, list) else [cluster_list],
        output_folder=f"{results_dirpath}/",
        outputfilename_prefix=cluster_header,
    )


def prefix_extension(filename, prefix):
    suffix = Path(filename).suffix
    return filename.replace(f"{suffix}", f"{prefix}{suffix}")


def main():
    """Parse command line arguments, then run NS-Forest function."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run NS-Forest functions.")
    parser.add_argument("h5ad_filepath", help="The dataset filepath")
    parser.add_argument(
        "-c",
        "--cluster-header",
        default="cell_type",
        help="The cluster header, default: 'cell_type'",
        metavar="CLUSTER HEADER",
    )
    parser.add_argument(
        "-l",
        "--cluster-list",
        default=[],
        help="The cluster list, default: all clusters",
        metavar="CLUSTER LIST",
    )
    parser.add_argument(
        "-d",
        "--results-dirpath",
        default=".",
        help="The directory path for results",
        metavar="DIR",
    )
    parser.add_argument(
        "-t",
        "--total-counts",
        default=0,
        type=int,
        help="Total counts after downsampling, default: 0, do not downsample",
        metavar="COUNTS",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--run-nsforest-with-preprocessing", action="store_true")
    group.add_argument("--preprocess-adata-file", action="store_true")
    group.add_argument("--downsample-adata-file", action="store_true")
    group.add_argument("--generate-scanpy-dendrogram", action="store_true")
    group.add_argument("--calculate-cluster-medians-per-gene", action="store_true")
    group.add_argument(
        "--calculate-binary-scores-per-gene-per-cluster", action="store_true"
    )
    group.add_argument("--run-nsforest-without-preprocessing", action="store_true")
    args = parser.parse_args()

    # Run the NS-Forest function
    if args.run_nsforest_with_preprocessing:
        run_nsforest_with_preprocessing(
            args.h5ad_filepath,
            cluster_header=args.cluster_header,
            cluster_list=args.cluster_list,
            results_dirpath=args.results_dirpath,
            total_counts=args.total_counts,
        )

    if args.preprocess_adata_file:
        preprocess_adata_file(
            args.h5ad_filepath,
            args.cluster_header,
            prefix_extension(args.h5ad_filepath, "_pp"),
        )

    if args.downsample_adata_file:
        downsample_adata_file(
            args.h5ad_filepath,
            args.total_counts,
            prefix_extension(args.h5ad_filepath, "_ds"),
        )

    if args.generate_scanpy_dendrogram:
        generate_scanpy_dendrogram(
            args.h5ad_filepath,
            args.cluster_header,
            prefix_extension(args.h5ad_filepath, "_gd"),
            results_dirpath=args.results_dirpath,
        )

    if args.calculate_cluster_medians_per_gene:
        calculate_cluster_medians_per_gene(
            args.h5ad_filepath,
            args.cluster_header,
            prefix_extension(args.h5ad_filepath, "_cm"),
        )

    if args.calculate_binary_scores_per_gene_per_cluster:
        calculate_binary_scores_per_gene_per_cluster(
            args.h5ad_filepath,
            args.cluster_header,
            prefix_extension(args.h5ad_filepath, "_cs"),
        )

    if args.run_nsforest_without_preprocessing:
        run_nsforest_without_preprocessing(
            args.h5ad_filepath,
            args.cluster_header,
            cluster_list=args.cluster_list,
            results_dirpath=args.results_dirpath,
        )


if __name__ == "__main__":
    main()
