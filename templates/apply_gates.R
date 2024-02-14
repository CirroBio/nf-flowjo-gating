#!/usr/bin/env Rscript

library(CytoML)
library(flowWorkspace)
library(ggcyto)
library(tidyr)
library(data.table)
library(matrixStats)
library(dplyr)
library(ggplot2)

ws <- open_flowjo_xml("${input_wsp}")
gh <- flowjo_to_gatingset(
  ws,
  includeGates=TRUE,
  name=1,
  execute=FALSE
)
fcs_files <- list.files(
  path=".",
  pattern="*.fcs",
  full.names=TRUE
)
cs <- load_cytoset_from_fcs(files=fcs_files)

print(cs)
gs <- gh_apply_to_cs(gh, cs)

nodelist <- gs_get_pop_paths(gs, path = "auto")
for(node in nodelist){
  if(node != "root"){
    autoplot(gs, node)
    ggsave(paste("gating.", node, ".pdf", sep=""))
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
  ggsave(paste(fix.pop.name(pop.name), prefix, "pdf", sep="."))
}

fix.pop.name <- function(pop.name){
  gsub("^_", "", gsub("/", "_", pop.name\$pop[1]))
}

write.summary.table <- function(dat, pop.name, prefix){
  write.table(
    dat,
    file=paste(
      fix.pop.name(pop.name),
      prefix,
      "csv",
      sep="."
    ),
    sep=",",
    quote=FALSE,
    row.names=FALSE
  )
}

# # Get the median fluorescence intensity for each marker
# group_map(
#   gs_pop_get_stats(gs, type=pop.MFI) %>% group_by(pop),
#   write.summary.MFI
# )

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
    paste("gating", metric, "csv", sep="."),
    row.names=FALSE,
    quote=FALSE
  )
}
