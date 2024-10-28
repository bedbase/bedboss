# R Script that generates statistics for a bedfile

library(GenomicDistributions)
library(GenomicDistributionsData)
library(GenomeInfoDb)
library(ensembldb)
library(optparse)
library(tools)
library(R.utils)
library(rjson)

trim <- IRanges::trim

option_list = list(
  make_option(c("--bedfilePath"), type="character", default=NULL, 
              help="full path to a BED file to process", metavar="character"),
  make_option(c("--fileId"), type="character", default=NULL,
              help="BED file ID to use for output files prefix", metavar="character"),
  make_option(c("--openSignalMatrix"), type="character",
              help="path to the open signal matrix required for the tissue specificity plot", metavar="character"),
  make_option(c("--digest"), type="character", default=NULL,
              help="digest of the BED file", metavar="character"),
  make_option(c("--outputFolder"), type="character", default="output",
              help="base output folder for results", metavar="character"),
  make_option(c("--genome"), type="character", default="hg38",
              help="genome reference to calculate against", metavar="character"),
  make_option(c("--ensdb"), type="character",
              help="path to the Ensembl annotation gtf file", metavar="character")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);    

if (is.null(opt$bedfilePath)) {
  print_help(opt_parser)
  stop("Bed file input missing.")
}

if (is.null(opt$fileId)) {
  print_help(opt_parser)
  stop("fileId input missing.")
}

if (is.null(opt$digest)) {
  print_help(opt_parser)
  stop("digest input missing.")
}

myPartitionList <- function(gtffile){
  features = c("gene", "exon", "three_prime_utr", "five_prime_utr")
  geneModels = getGeneModelsFromGTF(gtffile, features, TRUE)
  partitionList = genomePartitionList(geneModels$gene,
                                      geneModels$exon,
                                      geneModels$three_prime_utr,
                                      geneModels$five_prime_utr)
  
  return (partitionList) 
}


myChromSizes <- function(genome){
  if (requireNamespace(BSgm, quietly=TRUE)){
    library (BSgm, character.only = TRUE)
    BSG = eval(as.name(BSgm))
  } else {
    library (BSg, character.only = TRUE)
    BSG = eval(as.name(BSg))
  }
  chromSizesGenome = seqlengths(BSG)
  return(chromSizesGenome)
}

plotBoth <- function(plotId, g){
  pth = paste0(opt$outputFolder, "/", fileId, "_", plotId)
  print(paste0("Plotting: ", pth))
  ggplot2::ggsave(paste0(pth, ".png"), g, device="png", width=8, height=8, units="in")
  ggplot2::ggsave(paste0(pth, ".pdf"), g, device="pdf", width=8, height=8, units="in")
}

getPlotReportDF <- function(plotId, title){
  pth = paste0(opt$outputFolder, "/", fileId, "_", plotId)
  rel_pth = getRelativePath(pth, paste0(opt$outputFolder, "/../../../"))
  print(paste0("Writing plot json: ", rel_pth))
  newPlot = data.frame(
    "name"=plotId, 
    "title"=title, 
    "thumbnail_path"=paste0(rel_pth, ".png"), 
    "path"=paste0(rel_pth, ".pdf")
  )
  return(newPlot)
}


doItAall <- function(query, fileId, genome, cellMatrix) {
  plots = data.frame(stringsAsFactors=F)
  bsGenomeAvail = ifelse((requireNamespace(BSg, quietly=TRUE) | requireNamespace(BSgm, quietly=TRUE)), TRUE, FALSE)
  
  # check if json file exist for the input bed file
  meta_path = paste0(outfolder, "/", fileId, ".json")
  plot_path = paste0(outfolder, "/", fileId, "_plots.json")
  if (file.exists(meta_path)){
    bedmeta = fromJSON(file=meta_path)
    plots = fromJSON(file=plot_path) 
    plots = as.data.frame(do.call(rbind, plots))
  } else{
    plots = data.frame(stringsAsFactors=F)
  }

  run_plot = TRUE
  # TSS distance plot
  if (exists("bedmeta")){
    if ("median_TSS_dist" %in% names(bedmeta)){
      run_plot = FALSE
    } else {
      run_plot = TRUE
    }
  } else {
    run_plot = TRUE
  }
  if (run_plot){
    tryCatch(
      expr = {
        if (!(genome %in% c("hg19", "hg38", "mm10", "mm9")) && gtffile == "None"){
          message("Ensembl annotation gtf file not provided. Skipping TSS distance plot ... ")
        } else{
          if (genome %in% c("hg19", "hg38", "mm10", "mm9")) {
            TSSdist = calcFeatureDistRefTSS(query, genome)
            plotBoth("tssdist", plotFeatureDist( TSSdist, featureName="TSS"))
          } else {
            tss = getTssFromGTF(gtffile, TRUE)
            TSSdist = calcFeatureDist(query, tss)
            plotBoth("tssdist", plotFeatureDist( TSSdist, featureName="TSS"))
          }
        }
        if (exists("bedmeta")){
          tss <- list(median_TSS_dist = signif(median(abs(TSSdist), na.rm=TRUE), digits = 4))
          bedmeta = append(bedmeta, tss)
        }
        plots = rbind(plots, getPlotReportDF("tssdist", "Region-TSS distance distribution"))
        message("Successfully calculated and plot TSS distance.")
      },
      error = function(e){
        message('Caught an error in creating: TSS distance plot!')
        print(e)
      }
    )
  }
  
  
  
  # Chromosomes region distribution plot
  if (!exists("bedmeta") ){
    tryCatch(
      expr = {
        if (genome %in% c("mm39", "dm3", "dm6", "ce10", "ce11", "danRer10", "danRer10", "T2T")){
          chromSizes = myChromSizes(genome)
          genomeBins  = getGenomeBins(chromSizes)
          plotBoth("chrombins", plotChromBins(calcChromBins(query, genomeBins)))
        } else{
          plotBoth("chrombins", plotChromBins(calcChromBinsRef(query, genome)))
        }
        
        plots = rbind(plots, getPlotReportDF("chrombins", "Regions distribution over chromosomes"))
        message("Successfully calculated and plot chromosomes region distribution.")
      },
      error = function(e){
        message('Caught an error in creating: Chromosomes region distribution plot!')
        print(e)
      }
    ) 
  }
  
  
  
  
  # OPTIONAL: Plot GC content only if proper BSgenome package is installed. 
  if (exists("bedmeta")){
    if ("gc_content" %in% names(bedmeta)){
      run_plot = FALSE
    } else {
      run_plot = TRUE
    }
  } else {
      run_plot = TRUE
  }

  if (run_plot){
    if (bsGenomeAvail) {
      tryCatch(
        expr = {
          if (requireNamespace(BSgm, quietly=TRUE)){
            library (BSgm, character.only = TRUE)
            bsg = eval(as.name(BSgm))
            gcvec = calcGCContent(query, bsg)
          } else {
            library (BSg, character.only = TRUE)
            bsg = eval(as.name(BSg))
            gcvec = calcGCContent(query, bsg)
          }
          plotBoth("gccontent", plotGCContent(gcvec))
          if (exists("bedmeta")){
            gc_content <- list(gc_content = signif(mean(gcvec), digits = 4))
            bedmeta = append(bedmeta, gc_content)
          }
          plots = rbind(plots, getPlotReportDF("gccontent", "GC content"))
          message("Successfully calculated and plot GC content.")
        },
        error = function(e){
          message('Caught an error in creating: GC content plot!')
          print(e, gcvec)
        }
      ) 
    }
  }
  
  
  # Partition plots, default to percentages
  if (exists("bedmeta")){
    if ("exon_frequency" %in% names(bedmeta)){
      run_plot = FALSE
    } else {
      run_plot = TRUE
    }
  } else {
      run_plot = TRUE
  }  

  if (run_plot){
    tryCatch(
      expr = {
        if (!(genome %in% c("hg19", "hg38", "mm10")) && gtffile == "None"){
          message("Ensembl annotation gtf file not provided. Skipping partition plot ... ")
        } else {
          if (genome %in% c("hg19", "hg38", "mm10")) {
            gp = calcPartitionsRef(query, genome)
            plotBoth("partitions", plotPartitions(gp))
          } else {
            partitionList = myPartitionList(gtffile)
            gp = calcPartitions(query, partitionList)
            plotBoth("partitions", plotPartitions(gp))
          }
          plots = rbind(plots, getPlotReportDF("partitions", "Regions distribution over genomic partitions"))
          # flatten the result returned by the function above
          partiotionNames = as.vector(gp[,"partition"])
          partitionsList = list()
          for(i in seq_along(partiotionNames)){
            partitionsList[[paste0(partiotionNames[i], "_frequency")]] = 
              as.vector(gp[,"Freq"])[i]
            partitionsList[[paste0(partiotionNames[i], "_percentage")]] = 
              as.vector(gp[,"Freq"])[i]/length(query)	        
          }
          if (exists("bedmeta")){
            bedmeta = append(bedmeta, partitionsList)
          }
          message("Successfully calculated and plot regions distribution over genomic partitions.")
        }
      },
      error = function(e){
        message('Caught an error in creating: Partition plot!')
        print(e)
      }
    ) 
  }
 
  
  # Expected partition plots
  if (!exists("bedmeta") ){
    tryCatch(
      expr = {
        if (!(genome %in% c("hg19", "hg38", "mm10")) && gtffile == "None"){
          message("Ensembl annotation gtf file not provided. Skipping expected partition plot ... ")
        } else{
          if (genome %in% c("hg19", "hg38", "mm10")) {
            plotBoth("expected_partitions", plotExpectedPartitions(calcExpectedPartitionsRef(query, genome)))
          } else {
            partitionList = myPartitionList(gtffile)
            chromSizes = myChromSizes(genome)
            genomeSize = sum(chromSizes)
            plotBoth("expected_partitions", plotExpectedPartitions(calcExpectedPartitions(query, partitionList, genomeSize)))
          }
          plots = rbind(plots, getPlotReportDF("expected_partitions", "Expected distribution over genomic partitions"))
          message("Successfully calculated and plot expected distribution over genomic partitions.")
        }
      },
      error = function(e){
        message('Caught an error in creating: Expected partition plot!')
        print(e)
      }
    ) 
  }
 
  # Cumulative partition plots
  if (!exists("bedmeta") ){
    tryCatch(
      expr = {
        if (!(genome %in% c("hg19", "hg38", "mm10")) && gtffile == "None"){
          message("Ensembl annotation gtf file not provided. Skipping cumulative partition plot ... ")
        } else{
          if (genome %in% c("hg19", "hg38", "mm10")) {
            plotBoth("cumulative_partitions", plotCumulativePartitions(calcCumulativePartitionsRef(query, genome)))
          } else{
            partitionList = myPartitionList(gtffile)
            plotBoth("cumulative_partitions", plotCumulativePartitions(calcCumulativePartitions(query, partitionList)))
          }
          plots = rbind(plots, getPlotReportDF("cumulative_partitions", "Cumulative distribution over genomic partitions"))
          message("Successfully calculated and plot cumulative distribution over genomic partitions.")
        }
      },
      error = function(e){
        message('Caught an error in creating: Cumulative partition plot!')
        print(e)
      }
    ) 
  }
  
  # QThis plot
  if (exists("bedmeta")){
    if ("mean_region_width" %in% names(bedmeta)){
      run_plot = FALSE
    } else {
      run_plot = TRUE
    }
  } else {
      run_plot = TRUE
  }

  if (run_plot){
    tryCatch(
      expr = {
        widths = calcWidth(query)
        plotBoth("widths_histogram", plotQTHist(widths))
        if (exists("bedmeta")){
          mean_region_width <- list(mean_region_width = signif(mean(widths), digits = 4))
          bedmeta = append(bedmeta, mean_region_width)
        }
        plots = rbind(plots, getPlotReportDF("widths_histogram", "Quantile-trimmed histogram of widths"))
        message("Successfully calculated and plot quantile-trimmed histogram of widths.")
      },
      error = function(e){
        message('Caught an error in creating: Quantile-trimmed histogram of widths plot!')
        print(e, widths)
      }
    ) 
  }
  
  # Neighbor regions distance plots
  if (!exists("bedmeta") ){
    tryCatch(
      expr = {
        plotBoth("neighbor_distances", plotNeighborDist(calcNeighborDist(query)))
        plots = rbind(plots, getPlotReportDF("neighbor_distances", "Distance between neighbor regions"))
        message("Successfully calculated and plot distance between neighbor regions.")
      },
      error = function(e){
        message('Caught an error in creating: Distance between neighbor regions plot!')
        print(e)
      }
    ) 
  }
  
  
  # Tissue specificity plot if open signal matrix is provided
  if (!exists("bedmeta") ){
    if (cellMatrix == "None") {
      message("open signal matrix not provided. Skipping tissue specificity plot ... ")
    } else {
      tryCatch(
        expr = {
          plotBoth("open_chromatin", plotSummarySignal(calcSummarySignal(query, data.table::fread(cellMatrix))))
          plots = rbind(plots, getPlotReportDF("open_chromatin", "Cell specific enrichment for open chromatin"))
          message("Successfully calculated and plot cell specific enrichment for open chromatin.")
        },
        error = function(e){
          message('Caught an error in creating: Cell specific enrichment for open chromatin plot!')
          print(e)
        }
      )
    }
  }
 
  
  # Note: names of the list elements MUST match what's defined in: https://github.com/databio/bbconf/blob/master/bbconf/schemas/bedfiles_schema.yaml
  if (exists("bedmeta")){
    write(jsonlite::toJSON(bedmeta, pretty=TRUE), meta_path)
    write(jsonlite::toJSON(plots, pretty=TRUE, auto_unbox = TRUE), plot_path)
  } else {
    bedmeta = list(
      name=fileId,
      number_of_regions=length(query),
      mean_region_width=ifelse(exists('widths'), signif(mean(widths), digits = 4), NA),
      md5sum=opt$digest
    )
    if (exists('gcvec') && !isEmpty(gcvec)){
      gc_content <- list(gc_content = signif(mean(gcvec), digits = 4))
      bedmeta = append(bedmeta, gc_content)
    }
    if (exists('TSSdist') && !all(is.na(TSSdist))){
      tss <- list(median_TSS_dist = signif(median(abs(TSSdist), na.rm=TRUE), digits = 4))
      bedmeta = append(bedmeta, tss)
    }
    if (exists('partitionsList')){
      write(jsonlite::toJSON(c(bedmeta, partitionsList), pretty=TRUE), meta_path)
    } else {
      write(jsonlite::toJSON(c(bedmeta), pretty=TRUE), meta_path)
    }
    
    if (exists('plots')){
      write(jsonlite::toJSON(plots, pretty=TRUE, auto_unbox = TRUE), plot_path)
    }
  }
  }

# define values and output folder for doitall()
fileId = opt$fileId
bedPath = opt$bedfilePath
outfolder = opt$outputFolder
genome = opt$genome
cellMatrix = opt$openSignalMatrix
gtffile = opt$ensdb


# # build BSgenome package ID to check whether it's installed
# if ( startsWith(genome, "T2T")){
#   BSg = "BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0"
# } else {
#   if (startsWith(genome, "hg") | startsWith(genome, "grch")) {
#     orgName = "Hsapiens"
#   } else if (startsWith(genome, "mm") | startsWith(genome, "grcm")){
#     orgName = "Mmusculus"
#   } else if (startsWith(genome, "dm")){
#     orgName = "Dmelanogaster"
#   } else if (startsWith(genome, "ce")){
#     orgName = "Celegans"
#   } else if (startsWith(genome, "danRer")){
#     orgName = "Drerio"
#   }  else if (startsWith(genome, "TAIR")){
#     orgName = "Athaliana"
#   } else {
#     orgName = "Undefined"
#   }
#   BSg = paste0("BSgenome.", orgName , ".UCSC.", genome)
# }
#
# BSgm = paste0(BSg, ".masked")

# read bed file and run doitall()
query = LOLA::readBed(bedPath)
doItAall(query, fileId, genome, cellMatrix)
