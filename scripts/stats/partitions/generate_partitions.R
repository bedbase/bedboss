library(GenomicDistributions)
getGeneModels = function(refAssembly) {
  datasetId = paste0("geneModels_", refAssembly)

  if(refAssembly == "hg19"){

    geneModelsDataset = getReferenceData(refAssembly, tagline="geneModels_")

  } else if(refAssembly == "hg38"){

    if(!"GenomicDistributionsData" %in% utils::installed.packages()){
      stop(paste(datasetId, "not available in GenomicDistributions package",
                 "and GenomicDistributionsData package is not installed"))
    } else {
      geneModelsDataset = GenomicDistributionsData::geneModels_hg38()
    }

  } else if(refAssembly == "mm10"){

    if(!"GenomicDistributionsData" %in% utils::installed.packages()){
      stop(paste(datasetId, "not available in GenomicDistributions package",
                 "and GenomicDistributionsData package is not installed"))
    } else {
      geneModelsDataset = GenomicDistributionsData::geneModels_mm10()
    }

  } else if(refAssembly == "mm9"){

    if(!"GenomicDistributionsData" %in% utils::installed.packages()){
      stop(paste(datasetId, "not available in GenomicDistributions package",
                 "and GenomicDistributionsData package is not installed"))
    } else {
      geneModelsDataset = GenomicDistributionsData::geneModels_mm9()
    }

  } else {
    stop(paste(datasetId, "not available in GenomicDistributions package",
               "or GenomicDistributionsData package,",
               "please use getGeneModelsFromGTF() to get",
               "gene models."))
  }

  return(geneModelsDataset)
}

geneModels = getGeneModels("hg38")
partitionList = GenomicDistributions::genomePartitionList(geneModels$genesGR, 
                                          geneModels$exonsGR,
                                          geneModels$threeUTRGR, 
                                          geneModels$fiveUTRGR)

# partitionList$intron


calcPartitionsRef = function(query, refAssembly, bpProportion=FALSE){
    .validateInputs(list(query=c("GRanges", "GRangesList"),
                         refAssembly="character"))
    geneModels = getGeneModels(refAssembly)
    partitionList = genomePartitionList(geneModels$genesGR,
                                        geneModels$exonsGR,
                                        geneModels$threeUTRGR,
                                        geneModels$fiveUTRGR)
    message("Calculating overlaps...")
    return(calcPartitionsRef(query, partitionList, bpProportion=bpProportion))
}

path = "/home/bnt4me/virginia/repos/bedboss/scripts/stats/partitions/"

write.table(partitionList$intron, file = paste(path, "hg_38.intron.csv", sep=""), quote=F, row.names = F, col.names=FALSE, sep ='\t')
write.table(partitionList$promoterCore, file = paste(path, "hg_38.promoterCore.csv", sep=""), quote=F, row.names = F,sep="\t", col.names=FALSE)
write.table(partitionList$promoterProx, file = paste(path, "hg_38.promoterProx.csv", sep=""), quote=F, row.names = F,sep="\t", col.names=FALSE)
write.table(partitionList$exon, file = paste(path, "hg_38.exon.csv", sep=""), quote=F, row.names = F,sep="\t", col.names=FALSE)
write.table(partitionList$fiveUTR, file = paste(path, "hg_38.fiveUTR.csv", sep=""), quote=F, row.names = F,sep="\t", col.names=FALSE)
write.table(partitionList$threeUTR, file = paste(path, "hg_38.threeUTR.csv", sep=""), quote=F, row.names = F,sep="\t")
# 





