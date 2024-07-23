#!/usr/bin/env python3

import anndata as ad
import pandas as pd
from plotly.subplots import make_subplots
import plotly.express as px
from pathlib import Path
from scipy.stats import zscore


def read_adata() -> ad.AnnData:
    # Read in and concatenate any AnnData files with the name
    # input.*.h5ad
    return ad.concat(
        [
            ad.read_h5ad(f)
            for f in Path(".").glob("input.*.h5ad")
        ]
    )


def make_plot(cytometer: str, cytometer_serial: str, adata: ad.AnnData):
    # For each gated subset, make a plot of MFI and rCV
    # for each collector
    for subset in adata.obsm["count"].columns:

        # Make a plot with the number of cells in this subset
        # as a line plot
        fig = px.line(
            x=adata.obs["DATETIME"],
            y=adata.obsm["count"][subset],
            title=f"{subset} - Cell Count",
            labels={
                "x": "Collection Time",
                "y": "Cell Count"
            }
        )
        # Make the background transparent
        fig.update_layout(
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            hovermode="x"
        )
        # Write out to HTML
        html_output = Path(f"{cytometer}_{cytometer_serial}/{subset}/cell_count.html")
        html_output.parent.mkdir(parents=True, exist_ok=True)
        with open(html_output, "w") as f:
            f.write(fig.to_html(include_plotlyjs='cdn'))

        # Make sure that MFI and rCV data are available
        assert f"{subset}-MFI" in adata.layers
        assert f"{subset}-rCV" in adata.layers

        mfi = adata.to_df(f"{subset}-MFI").drop(columns=["Time"]).sort_index(axis=1)
        rcv = adata.to_df(f"{subset}-rCV").drop(columns=["Time"]).sort_index(axis=1)

        # Make a heatmap with the MFI and rCV data
        # The rows should be the collectors,
        # the columns should be the time points,
        # the colors should be z-score normalized values
        # (by row),
        # and the hover text should be the actual value

        fig = make_subplots(
            rows=2,
            cols=1,
            subplot_titles=("MFI", "rCV"),
            vertical_spacing=0.1,
            shared_xaxes=True,
            shared_yaxes=True
        )
        fig.add_heatmap(
            z=mfi.T.apply(zscore, axis=1).values,
            text=mfi.T.values,
            x=adata.obs["DATETIME"].apply(str).values,
            y=mfi.columns.values,
            row=1,
            col=1,
            colorscale="bluered",
            zmin=-2,
            zmax=2,
            hovertemplate="%{y} MFI (%{x}): %{text:.2f}",
            name="MFI",
            colorbar=dict(title="z-score")
        )
        fig.add_heatmap(
            z=rcv.T.apply(zscore, axis=1).values,
            text=rcv.T.values,
            x=adata.obs["DATETIME"].apply(str).values,
            y=rcv.columns.values,
            row=2,
            col=1,
            colorscale="bluered",
            zmin=-2,
            zmax=2,
            hovertemplate="%{y} rCV (%{x}): %{text:.2f}",
            name="rCV",
            colorbar=dict(title="z-score")
        )
        # Write out to HTML
        html_output = Path(f"{cytometer}_{cytometer_serial}/{subset}/heatmap.html")
        html_output.parent.mkdir(parents=True, exist_ok=True)
        with open(html_output, "w") as f:
            f.write(fig.to_html(include_plotlyjs='cdn'))

        for collector in mfi.columns:
            fig = make_subplots(specs=[[{"secondary_y": True}]])
            fig.add_scatter(
                x=adata.obs["DATETIME"],
                y=mfi[collector],
                mode="markers+lines",
                name="MFI",
                hovertemplate=None
            )
            fig.add_scatter(
                x=adata.obs["DATETIME"],
                y=rcv[collector],
                mode="markers+lines",
                name="rCV",
                secondary_y=True,
                hovertemplate=None
            )
            # Make the background transparent
            fig.update_layout(
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                hovermode="x",
                title_text=collector,
                xaxis_title="Collection Time",
                yaxis_title="MFI",
            )
            fig.update_yaxes(title_text="rCV", secondary_y=True)
            # Set up the output folder
            html_output = Path(f"{cytometer}_{cytometer_serial}/{subset}/channels/{collector}.html")
            html_output.parent.mkdir(parents=True, exist_ok=True)
            # Write out the HTML file with the figure
            with open(html_output, "w") as f:
                f.write(fig.to_html(include_plotlyjs='cdn'))


def main():
    adata = read_adata()

    # Parse the DATE column as a datetime with the format DD-Month-YYYY
    # while also adding in the data from the BTIM column
    # which is of the format "HH:MM:SS"
    adata.obs["DATETIME"] = pd.to_datetime(
        adata.obs.apply(
            lambda r: r["DATE"] + " " + r["BTIM"],
            axis=1
        ),
        format="%d-%b-%Y %H:%M:%S"
    )

    adata = adata[adata.obs["DATETIME"].sort_values().index]

    # For each instrument, make a set of plots
    for (cytometer, cytometer_serial), subset in adata.obs.groupby(["CYT", "CYTNUM"]):
        make_plot(cytometer, cytometer_serial, adata[subset.index])


if __name__ == '__main__':
    main()
