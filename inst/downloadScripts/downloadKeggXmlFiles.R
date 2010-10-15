library(R.utils)
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

rootPath <- "ftp.genome.jp/pub/kegg/xml/kgml"
organism <- "hsa"

localPath <- "networkData"
localPath <- Arguments$getWritablePath(localPath)

if (!exists("force", mode="character")) {
  force <- FALSE
}

for (tag in c("non-metabolic", "metabolic")) {
  path <- file.path(rootPath, tag, "organisms", organism)
  srcPath <- file.path("ftp:/", path, "")
  tgtPath <- file.path(localPath, path)
  tgtPath <- Arguments$getWritablePath(tgtPath)

  verbose && cat(verbose, "Remote directory:")
  verbose && print(verbose, srcPath)

  verbose && enter(verbose, "Writing directory contents to INDEX file")
  tgtPathname <- file.path(tgtPath, "INDEX")
  downloadFile(srcPath, tgtPathname, method="wget", verbose=verbose)
  verbose && exit(verbose)
  
  fl <- readLines(tgtPathname)
  rm(tgtPathname)
  gg <- grep("File", fl)
  dat <- fl[gg]
  pattern <-   ".*<a href=\\\"(.*)\\\">(.*\\.xml).*"
  srcPathnames <- gsub(pattern, "\\1", dat)
  tgtFilenames <- gsub(pattern, "\\2", dat)
  verbose && cat(verbose, "Contents:")
  verbose && str(verbose, tgtFilenames)

  for (kk in seq(along=srcPathnames)) {
    srcPathname <- srcPathnames[kk]
    tgtFilename <- tgtFilenames[kk]
    tgtPathname <- file.path(tgtPath, tgtFilename)

    downloadFile(srcPathname, tgtPathname, skip=!force, overwrite=force, verbose=verbose)
    ## TODO: remove local files that are not present on the remote server?
  }
}

############################################################################
# HISTORY:
# 2010-09-14
# o Created.
############################################################################
