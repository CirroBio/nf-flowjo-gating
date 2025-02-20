#!/usr/bin/env Rscript

library(CytoML)
library(flowWorkspace)
library(ggcyto)
library(tidyr)
library(data.table)
library(matrixStats)
library(dplyr)
library(ggplot2)

# Set up a function to create any parent folders needed for a particular path
make_parent_folder <- function(fp){
  items <- strsplit(fp, "/")[[1]]
  if(length(items) > 1){
    folder <- paste(items[1:(length(items)-1)], collapse="/")
    print(folder)
    if(!file.exists(folder)){
      print(paste("Creating parent folder", folder))
      dir.create(folder, recursive=TRUE)
    }
  }
}

print("Creating directory gating/")
dir.create("gating", showWarnings = FALSE)

print("Opening gating file ${input_wsp}")
ws <- open_flowjo_xml("${input_wsp}")

print("Converting to gating set")
gh <- flowjo_to_gatingset(
  ws,
  includeGates=TRUE,
  name=1,
  execute=FALSE
)

print("Listing FCS input files (*.fcs)")
fcs_files <- list.files(
  path=".",
  pattern="*.fcs",
  full.names=TRUE
)

print("Writing out metadata.csv")
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

print("Loading cytoset from FCS files")
cs <- load_cytoset_from_fcs(files=fcs_files)

print("Applying gates to input data")
gs <- gh_apply_to_cs(gh, cs)

nodelist <- gs_get_pop_paths(gs, path = "auto")
for(node in nodelist){
  if(node != "root"){
    print(paste("Plotting node", node))
    autoplot(gs, node)
    output_fp = paste("gating/", node, ".pdf", sep="")
    make_parent_folder(output_fp)
    print(paste("Writing out", output_fp))
    ggsave(output_fp)
  }
}

pop.MFI <- function(fr){
  pd <- pData(parameters(fr))
  pd <- data.table(pd)
  chnls  <- pd[, name]

  res <- colMedians(exprs(fr)[, chnls, drop = FALSE])
  res
}


pop.rCV <- function(fr){
  pd <- pData(parameters(fr))
  pd <- data.table(pd)
  chnls  <- pd[, name]

  res <- apply(exprs(fr)[, chnls, drop = FALSE], 2, quantile, probs= c(0.25, 0.5, 0.75))
  res <- (res["75%",] - res["25%",]) / res["50%",]
  res
}

write.summary.MFI <- function(dat, pop.name){
  write.summary(dat, pop.name, "MFI", zscore=TRUE)
}
write.summary.rCV <- function(dat, pop.name){
  write.summary(dat, pop.name, "rCV", zscore=FALSE)
}

write.summary <- function(dat, pop.name, prefix, zscore=TRUE){
  write.summary.table(dat, pop.name, prefix)
  write.summary.heatmap(dat, pop.name, prefix, zscore=zscore)
}

write.summary.heatmap <- function(dat, pop.name, prefix, zscore=TRUE){
  plot_df <- pivot_longer(
    dat,
    colnames(dat)[2:ncol(dat)]) %>% filter(name != "Time")
  
  setDT(plot_df)[, Z.Score := scale(value), name]
  
  g = ggplot(plot_df,aes(x=sample, y=name)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab("") +
    ylab("Channel")
  
  if(zscore){
    g + 
    geom_tile(aes(fill=Z.Score)) +
    ggtitle(paste(prefix, "(z-score)", pop.name))

  }else{
    g + 
    geom_tile(aes(fill=value)) +
    ggtitle(paste(prefix, pop.name))
  }
  output_folder = paste("subset_summaries", fix.pop.name(pop.name), sep="/")
  dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
  output_fp = paste(output_folder, paste(prefix, "pdf", sep="."), sep="/")
  print(paste("Writing out", output_fp))
  make_parent_folder(output_fp)
  ggsave(output_fp)
}

fix.pop.name <- function(pop.name){
  gsub("^_", "", gsub("/", "_", pop.name\$pop[1]))
}

write.summary.table <- function(dat, pop.name, prefix){
  output_folder = paste("subset_summaries", fix.pop.name(pop.name), sep="/")
  dir.create(output_folder, showWarnings = FALSE, recursive=TRUE)
  output_fp = paste(
    output_folder,
    paste(
      prefix,
      "csv",
      sep="."
    ),
    sep="/"
  )
  print(paste("Writing out", output_fp))
  make_parent_folder(output_fp)
  write.table(
    dat,
    file=output_fp,
    sep=",",
    quote=FALSE,
    row.names=FALSE
  )
}

# Get the median fluorescence intensity for each marker
group_map(
  gs_pop_get_stats(gs, type=pop.MFI) %>% group_by(pop),
  write.summary.MFI
)

# Get the robust coefficient of variation for each marker
gs_pop_get_stats(gs, type=pop.rCV)
group_map(
  gs_pop_get_stats(gs, type=pop.rCV) %>% group_by(pop),
  write.summary.rCV
)

for(metric in c("percent", "count")){
  df <- gs_pop_get_stats(gs, type=metric)
  wide <- pivot_wider(
    df,
    id_cols=c("sample"),
    names_from="pop",
    values_from=all_of(metric)
  )
  write.csv(
    wide,
    paste("gating", paste(metric, "csv", sep="."), sep="/"),
    row.names=FALSE,
    quote=FALSE
  )
}
