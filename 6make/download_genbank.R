acc=paste("AB115",403:404,sep="")

download_genbank <- function(accession) {
  for (acc in accession) {
    URL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                 paste(acc, collapse = ","), "&rettype=gb&retmode=text",
                 sep = "")
    print(URL)
    #utils::download.file(url = URL,mode='libcurl', destfile = paste0(acc, ".gb"), quiet = TRUE)
    
    cmd = paste('curl', paste0("\'", URL, "\'"), '-o', paste0(acc, ".gb"))
    print(cmd)
    system(cmd)
  }
}

download_genbank(acc)