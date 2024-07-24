import argparse
import os

import scanpy as sc

import nsforest as ns
from nsforest import nsforesting


def run_nsforest_on_file(h5ad_filepath, cluster_header="cell_type", total_counts=5000):
    """Run NSForest using the specified dataset filepath, and cluster_header.

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


if __name__ == "__main__":
    main()
