#!/usr/bin/env Rscript

library(CytoML)
library(flowWorkspace)
library(ggcyto)
install.packages("tidyr")
library(tidyr)

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
    ggsave(paste(node, ".pdf", sep=""))
  }
}

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
    paste(metric, "csv", sep="."),
    row.names=FALSE,
    quote=FALSE
  )
}
