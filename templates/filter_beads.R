#!/usr/bin/env Rscript

library(flowCore)
library(ggplot2)

# Output files that will be created:
#   metadata.csv - Details on how the information in each file was collected
#   subset_summaries/root/MFI.csv - MFI for each channel, before filtering
#   subset_summaries/root/rCV.csv - rCV for each channel, before filtering
#   subset_summaries/beads/MFI.csv - MFI for each channel, after filtering
#   subset_summaries/beads/rCV.csv - rCV for each channel, after filtering
#   gating/count.csv - Number of cells - columns are root and beads
#   gating/percent.csv - Percent of total cells - columns are root and beads

dir.create("gating", showWarnings = FALSE)
dir.create("subset_summaries", showWarnings = FALSE)
dir.create("subset_summaries/root", showWarnings = FALSE)
dir.create("subset_summaries/beads", showWarnings = FALSE)

fcs_files <- list.files(
    path=".",
    pattern="*.fcs",
    full.names=TRUE
)

# When adding the filename to the output, remove the ./ prefix
fix_filename <- function(fp){
    if(substr(fp, 1, 2) == "./"){
        return(substr(fp, 3, nchar(fp)))
    }
}

# Write out the metadata information for all files
write.csv(
    do.call(
        rbind,
        lapply(
            fcs_files,
            function(fp){
                keyword(read.FCS(fp))[
                    c(
                        "\$FIL",
                        "\$TOT",
                        "CREATOR",
                        "TUBE NAME",
                        "\$SRC",
                        "EXPERIMENT NAME",
                        "GUID",
                        "\$DATE",
                        "\$BTIM",
                        "\$ETIM",
                        "\$CYT",
                        "CYTNUM",
                        "EXPORT TIME"
                    )
                ]
            }
        )
    ),
    "metadata.csv",
    quote=FALSE,
    row.names=FALSE
)

# Set up the data structures which will ultimately be written out
root.MFI <- data.frame()
root.rCV <- data.frame()
beads.MFI <- data.frame()
beads.rCV <- data.frame()
count <- data.frame()
percent <- data.frame()

# Set the gates using a multiplier on the interquartile range
gate_size_multiplier <- 4
print(paste("Gate size:", gate_size_multiplier, "times the IQR"))

# Iterate over each FCS file
for(fp in fcs_files){

    # Read in the FCS file
    print(paste("Reading", fp))
    fcs <- read.FCS(fp)

    # Build the gate for beads
    gates <- list()
    gates[["filterId"]] <- "beads"
    pdf(paste("gating/", fp, ".beads.pdf", sep=""))
    for(channel in colnames(fcs)){
        if(channel == "Time") next
        # Get the values for all beads for this channel
        chVals <- fcs@exprs[,channel]
        # Calculate the median and quartiles
        sumVals <- summary(chVals)
        print(channel)
        print(sumVals)
        # Add a gate around the median which triples the interquartile range
        lowerIQR <- sumVals[["Median"]] - sumVals[["1st Qu."]]
        upperIQR <- sumVals[["3rd Qu."]] - sumVals[["Median"]]
        gates[[channel]] <- c(sumVals[["Median"]] - (gate_size_multiplier * lowerIQR), sumVals[["Median"]] + (gate_size_multiplier * upperIQR))
        print(paste("Keeping", gates[[channel]][1], "to", gates[[channel]][2]))
        # Make a plot
        plot(
          ggplot(
            data.frame(vals=chVals),
            aes(x=vals)
          )
          + geom_histogram(bins=150)
          + ggtitle(channel)
          + xlim(gates[[channel]][1], gates[[channel]][2])
        )
    }
    dev.off()
    # Build a gate
    beadsGate <- do.call(rectangleGate, gates)
    # Apply the gates, filtering the input data
    beads <- flowCore::filter(fcs, beadsGate)
    print("Filtering Complete")
    print(summary(beads))
    
    # Save the count and percent of cells
    count <- rbind(
        count,
        data.frame(
        file=fix_filename(fp),
        root=nrow(fcs@exprs),
        beads=sum(beads@subSet)
        )
    )
    percent <- rbind(
        percent,
        data.frame(
        file=fix_filename(fp),
        root=100,
        beads=100 * sum(beads@subSet) / nrow(fcs@exprs)
        )
    )

    # Calculate the robust coefficient of variation
    calc_rCV <- function(vals){
        summ <- summary(vals)
        (summ[["3rd Qu."]] - summ[["1st Qu."]]) / summ[["Median"]]
    }
    
    # Compute a summary for each channel using only
    # those measurements which pass the filter
    root.MFI_vals <- data.frame(file=fix_filename(fp))
    root.rCV_vals <- data.frame(file=fix_filename(fp))
    beads.MFI_vals <- data.frame(file=fix_filename(fp))
    beads.rCV_vals <- data.frame(file=fix_filename(fp))
    for(channel in colnames(fcs)){
        if(channel == "Time") next
        root.MFI_vals[[channel]] <- mean(fcs@exprs[, channel])
        beads.MFI_vals[[channel]] <- mean(fcs@exprs[beads@subSet, channel])
        root.rCV_vals[[channel]] <- calc_rCV(fcs@exprs[, channel])
        beads.rCV_vals[[channel]] <- calc_rCV(fcs@exprs[beads@subSet, channel])
    }
    # Add to the running total for all files
    root.MFI <- rbind(root.MFI, root.MFI_vals)
    root.rCV <- rbind(root.rCV, root.rCV_vals)
    beads.MFI <- rbind(beads.MFI, beads.MFI_vals)
    beads.rCV <- rbind(beads.rCV, beads.rCV_vals)

}

# Write out all of the combined results
write.csv(count, "gating/count.csv", row.names = FALSE, quote = FALSE)
write.csv(percent, "gating/percent.csv", row.names = FALSE, quote = FALSE)
write.csv(root.MFI, "subset_summaries/root/MFI.csv", row.names = FALSE, quote = FALSE)
write.csv(root.rCV, "subset_summaries/root/rCV.csv", row.names = FALSE, quote = FALSE)
write.csv(beads.MFI, "subset_summaries/beads/MFI.csv", row.names = FALSE, quote = FALSE)
write.csv(beads.rCV, "subset_summaries/beads/rCV.csv", row.names = FALSE, quote = FALSE)
