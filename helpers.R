# Collect all genes related to certain cellular function in GO.db using a regex. GO.db has to be loaded! This is species independent
grepGo <- function(term){
  goTerms <- unlist(Term(GOTERM))
  grepTerms <- grep(term, goTerms)
  goTerms[grepTerms]
} 

# Use output of grepGO to collect all genes related to certain grepped term or other GO ID list named with the GO ID. Returns the unique genes attached to the terms, not to the children terms.
collectGenesGrepGo_mouse <- function(grepGo, ...){
  select(org.Mm.eg.db, keys = names(grepGo), keytype = "GO", columns = c("GO", "SYMBOL", "ENTREZID"))
}

collectGenesGrepGo_human <- function(grepGo, ...){
  select(org.Hs.eg.db, keys = names(grepGo), keytype = "GO", columns = c("GO", "SYMBOL", "ENTREZID"))
}

# Combine the grep of GO terms and collection of genes in a single function. Return a vector of unique gene symbols or entrez ids related to the GO term.
collectSymbolGoTerm_mouse <- function(term){
  tmpGrepGo <- grepGo(term)
  if(length(tmpGrepGo>0)){
  tmpGenes <- collectGenesGrepGo_mouse(tmpGrepGo)
  } else {
    stop("This query did not return GO terms!")
  }
    return(unique(as.vector(tmpGenes[,"SYMBOL"])))
}

collectEntrezGoTerm_mouse <- function(term){
  tmpGrepGo <- grepGo(term)
  if(length(tmpGrepGo>0)){
    tmpGenes <- collectGenesGrepGo_mouse(tmpGrepGo)
  } else {
    stop("This query did not return GO terms!")
  }
    return(unique(as.vector(tmpGenes[,"ENTREZID"])))
}

# Map entrez to symbol, this can return multiple gene symbols per entrez id but this function only returns the unique names!
convertEntrez <- function(x){
  tmp <- select(org.Mm.eg.db, keys = x, keytype = "ENTREZID", columns =c("ENTREZID", "SYMBOL"))
  unique(as.vector(tmp[,"SYMBOL"]))
}
