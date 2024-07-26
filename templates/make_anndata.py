#!/usr/bin/env python3

import anndata as ad
import pandas as pd
from pathlib import Path


def main():

    # Read in the metadata file
    obs = (
        pd.read_csv("metadata.csv", index_col=0)
        .rename(columns=lambda cname: cname.strip('\$'))
    )

    # Read in all of the data as layers
    layers = {
        f"{folder.name}-{dtype}": (
            pd.read_csv(
                folder / f"{dtype}.csv",
                index_col=0
            )
            # Align the data with the metadata
            .reindex(index=obs.index)
        )
        for folder in Path("subset_summaries").iterdir()
        if folder.is_dir() and folder.name != "gating"
        for dtype in ["MFI", "rCV"]
    }

    # Read in the number and percentage of cells
    # in each gate as obsm
    obsm = {
        dtype: (
            pd.read_csv(
                f"gating/{dtype}.csv",
                index_col=0
            )
            .rename(
                columns=lambda cname: (
                    cname
                    .strip('/')
                    .replace('/', '_')
                )
            )
            # Align the data with the metadata
            .reindex(index=obs.index)
        )
        for dtype in ["count", "percent"]
    }

    # Validate the structure. Each subset which has
    # a set of counts and percentages should have
    # measurements for rCV and MFI as well.
    for subset in obsm["count"].columns.values:
        for dtype in ["MFI", "rCV"]:
            if f"{subset}-{dtype}" not in layers:
                raise ValueError(
                    f"Missing {dtype} data for subset {subset}"
                )

    # Make an AnnData object
    adata = ad.AnnData(
        X=layers["root-MFI"],
        layers=layers,
        obs=obs,
        obsm=obsm
    )

    adata.write_h5ad("output.h5ad")


if __name__ == '__main__':
    main()
