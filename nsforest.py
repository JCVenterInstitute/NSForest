#!/usr/bin/env python

import argparse
import os

import scanpy as sc

import nsforest as ns
from nsforest import nsforesting


def run_nsforest_on_file(h5ad_filepath, cluster_header="cell_type", total_counts=5000):
    """Run NSForest using the specified dataset filepath, and
    cluster_header.

    Parameters
    ----------
    h5ad_filepath : str
        The dataset filepath
    cluster_header : str
        The cluster header
    total_counts : int
        Total counts after downsampling

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
    # TODO: Improve
    results_dirpath = "."
    pp_h5ad_filename = f"pp_{h5ad_filename}"
    pp_h5ad_filepath = f"{results_dirpath}/{pp_h5ad_filename}"

    # Run NSForest if results do not exist
    if not os.path.exists(pp_h5ad_filepath):

        print(f"Loading unprocessed AnnData file: {h5ad_filename}")
        up_adata = sc.read_h5ad(h5ad_filepath)

        print("Calculating QC metrics")
        up_metrics = sc.pp.calculate_qc_metrics(up_adata)
        if up_metrics[1]["total_counts"].sum() > total_counts:
            print("Downsampling unprocessed AnnData file")
            ds_adata = sc.pp.downsample_counts(
                up_adata, total_counts=total_counts, copy=True
            )
        else:
            ds_adata = up_adata  # No need to copy

        print("Generating scanpy dendrogram")
        # Dendrogram order is stored in
        # `pp_adata.uns["dendrogram_cluster"]["categories_ordered"]`
        pp_adata = up_adata.copy()
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
            output_folder=f"{results_dirpath}/",
            outputfilename_prefix=cluster_header,
        )

    else:
        print(f"Completed NSForest for preprocessed AnnData file: {pp_h5ad_filename}")


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
        return out_adata

    else:
        return inp_adata

    print(f"Saving output AnnData file: {out_adata_file}")
    out_adata.write_h5ad(out_adata_file)


def generate_scanpy_dendrogram(
    inp_adata_file, cluster_header, out_dendrogram_dir, out_adata_file
):
    """Plots a dendrogram of the categories defined in grouping by
    cluster_header.

    Parameters
    ----------
    inp_adata_file : str
        Input AnnData file name
    cluster_header : str
        Column in `adata.obs` storing cell annotation
    out_dendrogram_dir : str
        Output directory, created if doesn't exist
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
    #   `pp_adata.uns["dendrogram_cluster"]["categories_ordered"]`
    out_adata = inp_adata.copy()
    out_adata.obs[cluster_header] = out_adata.obs[cluster_header].astype(str)
    out_adata.obs[cluster_header] = out_adata.obs[cluster_header].astype("category")
    out_adata = ns.pp.dendrogram(
        out_adata,
        cluster_header,
        save=False,
        output_folder=out_dendrogram_dir,
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


def run_nsforest(inp_adata_file, cluster_header, out_csv_dir):
    """Performs the main NS-Forest algorithm to find a list of
    NS-Forest markers for each `cluster_header`.

    Parameters
    ----------
    inp_adata_file : str
        Input AnnData file name
    cluster_header : str
        Column in `adata.obs` storing cell annotation
    out_csv_dir : str
        Output CSV directory

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
        output_folder=f"{out_csv_dir}/",
        outputfilename_prefix=cluster_header,
    )


def main():
    """Parse command line arguments, then run NS-Forest."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run NS-Forest.")
    parser.add_argument("h5ad_filepath", help="The dataset filepath")
    parser.add_argument(
        "-c",
        "--cluster-header",
        help="The cluster header, default: 'cell_type'",
        metavar="CLUSTER",
    )
    parser.add_argument(
        "-t",
        "--total-counts",
        type=int,
        help="Total counts after downsampling, default: 5000",
        metavar="COUNTS",
    )
    args = parser.parse_args()

    # Run NS-Forest
    if args.cluster_header or args.total_counts:
        if not args.cluster_header:
            run_nsforest_on_file(args.h5ad_filepath, args.total_counts)

        elif not args.total_counts:
            run_nsforest_on_file(args.h5ad_filepath, args.cluster_header)

        else:
            run_nsforest_on_file(
                args.h5ad_filepath, args.cluster_header, args.total_counts
            )

    else:
        run_nsforest_on_file(args.h5ad_filepath)

        inp_adata_file = args.h5ad_filepath
        total_counts = 5000
        out_adata_file = args.h5ad_filepath.lower().replace(".h5ad", "_1.h5ad")

        downsample_adata_file(inp_adata_file, total_counts, out_adata_file)

        inp_adata_file = out_adata_file
        cluster_header = "cell_type"
        out_dendrogram_dir = ""
        out_adata_file = inp_adata_file.replace("_1.", "_2.")

        generate_scanpy_dendrogram(
            inp_adata_file, cluster_header, out_dendrogram_dir, out_adata_file
        )

        inp_adata_file = out_adata_file
        cluster_header = "cell_type"
        out_adata_file = inp_adata_file.replace("_2.", "_3.")

        calculate_cluster_medians_per_gene(
            inp_adata_file, cluster_header, out_adata_file
        )

        inp_adata_file = out_adata_file
        cluster_header = "cell_type"
        out_adata_file = inp_adata_file.replace("_3.", "_4.")

        calculate_binary_scores_per_gene_per_cluster(
            inp_adata_file, cluster_header, out_adata_file
        )

        inp_adata_file = out_adata_file
        cluster_header = "cell_type"
        out_csv_dir = "."

        run_nsforest(inp_adata_file, cluster_header, out_csv_dir)


if __name__ == "__main__":
    main()
