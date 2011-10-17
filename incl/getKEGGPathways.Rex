library("Rgraphviz")
library("KEGGgraph")

## example of KGML files
path <- system.file("extdata", package="KEGGgraph")
grList <- getKEGGPathways(path=path, verbose=TRUE)
print(grList)

graph <- grList[[1]]
plotKEGGgraph(graph)

\dontrun{
  ## Download all human KGML pathways locally
  pathname <- system.file("downloadScripts", "downloadKeggXmlFiles.R", package="DEGraph")
  source(pathname)

  ## Load some of them
  grList <- getKEGGPathways(pattern="040", verbose=TRUE)
  print(grList)

  graph <- grList[[1]]
  plotKEGGgraph(graph)
}
