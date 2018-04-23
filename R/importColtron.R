#' read_coltron_cliques
#' @export
read_coltron_cliques <- function(path) {
  clique_file <- Sys.glob(paste0(path,"/*CLIQUES_ALL.txt"))
  if (!file.exists(clique_file)) stop("Clique file is missing.")
  clique_file %>%
    readLines %>%
    stringr::str_trim() %>%
    lapply(stringr::str_split_fixed,n=Inf,pattern='\t') %>%
    lapply(as.vector) %>%
    list_to_cliquelist()
}


#' read_coltron_subpeaks
#' @export
read_coltron_subpeaks <- function(path) {
  subpeak_file <- Sys.glob(paste0(path,"/*subpeaks.bed"))
  if (!file.exists(subpeak_file)) stop("Subpeak file is missing.")
  subpeak_file %>% rtracklayer::import()
}


#' read_coltron_tfbs
#' @export
read_coltron_tfbs <- function(path) {
  motif_files <- Sys.glob(paste0(path,"/motifBED/*motifs.bed"))
  if (!all(file.exists(motif_files))) stop("Motif .bed paths are weird.")
  motif_names <-  motif_files %>%
    stringr::str_split('/') %>%
    lapply(tail,1) %>%
    stringr::str_split("_") %>%
    lapply(head,1) %>%
    unlist
  names(motif_files) <- motif_names
  motif_files %>% lapply(rtracklayer::import) %>%
    GRangesList
}

#' import_coltron_sample
#' @export
import_coltron_sample <- function(path, name, bam=NULL) {
  if (!file.exists(path)) stop("This COLTRON directory doesn't exist.")
  message(paste0("Importing ", name, "... "),appendLF = F)
  tictoc::tic()
  subpeaks <- read_coltron_subpeaks(path)
  cliques <- read_coltron_cliques(path)
  tfbs <- read_coltron_tfbs(path)

  crcv <- CRCView(subpeaks, cliques, prior_tfbs=tfbs,
                  sample=name, bampath=bam)
  tictoc::toc()
  return(crcv)
}

#' create_coltron_experiment
#'
#' supply a data.frame containing at least:
#' SAMPLE, CONDITION, COLTRONDIR, BAM
#'
#' QUANTSITES is an optional field for specifying atac regions to
#' search for motifs in and quantify signal. Defaults to coltron's subpeaks
#' bed file otherwise.
#' QUANTMODE is an optional field for specifying whether to count reads like an atac
#' or a chipseq experiment. 'ATAC' or 'READS'.
#'
#' Result inherits from RangedSummarizedExperiment and is
#' compatible with regular chromvar workflow.
#'
#' @export
create_coltron_experiment <- function(metadata, quantsites="SUBPEAKS",
                                      resizeWidth = 1000, quantmode="READS",
                                      nthreads = 1, cohort_name="Coltron",
                                      genome='BSgenome.Hsapiens.UCSC.hg19') {
  stopifnot(class(metadata) == "data.frame")
  fields <- c("SAMPLE", "CONDITION", "COLTRONDIR", "BAM")
  stopifnot(fields %in% colnames(metadata))
  metadata %>%
    set_rownames(metadata$SAMPLE) -> metadata
  if (!("QUANTSITES" %in% colnames(metadata))) metadata$QUANTSITES <- quantsites
  metadata$BAM %>% as.vector %>% file.exists() %>% all %>% stopifnot()
  vanilla_crcview_list <- list()
  for(i in 1:nrow(metadata)) {
    dir <- as.vector(metadata[i,"COLTRONDIR"])
    nm <- as.vector(metadata[i,"SAMPLE"])
    bam <- as.vector(metadata[i,"BAM"])
    qs <- as.vector(metadata[i,"QUANTSITES"])
    import_coltron_sample(dir,nm,bam) -> crcv
    if (qs != "SUBPEAKS") {
      if(!file.exists(qs)) stop("QUANTSITE: Use existing .bed or default.")
      qs_gr <- rtracklayer::import(qs)
      CRCView(qs_gr, crcv@cliques, prior_tfbs=crcv@tfbs,
                          sample=crcv@name, bampath=crcv@bam) -> crcv
    }
    vanilla_crcview_list <- list(vanilla_crcview_list,crcv)
  }
  crcvlist <- CRCViewList(vanilla_crcview_list)

  all_candidate_sites <- quantsites(crcvlist) %>%
    unlist %>%
    resize(width = resizeWidth, fix="center") %>%
    set_names(.,as.character(.)) %>%
    unique %>% sort
  bams <- bam(crcvlist)
  message("Quantifying signal...", appendLF = F)
  tictoc::tic()
  capture.output(
  if (quantmode=="ATAC") {
    quantifyCutsites(all_candidate_sites,
                     bamlist = bams, nthreads =  nthreads)  -> cts
  } else {
    quantifyReads(all_candidate_sites,
                     bamlist = bams, nthreads =  nthreads)  -> cts
  },
  file="/dev/null")
  tictoc::toc()
  colnames(cts) <- names(bams)
  metadata$depth <- colSums(cts)
  rse <- SummarizedExperiment(assays = list(counts=cts),
                              rowRanges = all_candidate_sites,
                              colData = metadata[c("CONDITION","depth")])
  metadata(rse) <- metadata
  message("Finding non-overlapping peaks... ",appendLF = F)
  tictoc::tic()
  rse %<>% GenomeInfoDb::sortSeqlevels() %>% sort %>%
    chromVAR::filterPeaks(non_overlapping = T)
  tictoc::toc()
  message("Adding GC bias... ",appendLF = F)
  tictoc::tic()
  rse %<>% chromVAR::addGCBias(genome=genome)
  tictoc::toc()
  CRCExperiment(rse, crclist = crcvlist, cohort = cohort_name)
}




