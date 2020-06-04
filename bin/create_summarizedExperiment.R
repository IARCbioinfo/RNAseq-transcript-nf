args = commandArgs(trailingOnly=TRUE)
annot=args[1]
provider=args[2]
version=args[3]
genome=args[4]
organism=args[5]
gtf.path=args[6]
suffix=args[7]
fasta=args[8]

# print arguments
print(args)
#
samples = list.dirs(".")
samples = gsub("./", "", samples[samples!="."])

# transcriptome R file creation
library(tximeta)
library(tximport)
library(readr)
library(dplyr)
library(SummarizedExperiment)

# list files and transcript/gene link
STfiles = list.files(".",pattern = "t_data.ctab",recursive = T,full.names = T)
STnames = list.dirs(".")
STnames = gsub("./", "", STnames[STnames!="."])
coldata <- data.frame(files=STfiles, names=STnames, stringsAsFactors=FALSE)
tx2gene = read.table(STfiles[1],h=T)[,c("t_name","gene_id")]
colnames(tx2gene) = c("TXNAME","GENEID")

# import files and create summarizedExperiment object
txim = tximeta(coldata = coldata,type = "stringtie",txOut = T,tx2gene=tx2gene)
txim.genes = tximeta(coldata = coldata,type = "stringtie",txOut = F,tx2gene=tx2gene)

# create Txome annotation for our annot
gr <- rtracklayer::import(annot) #, format = "GFF", colnames = colnames, 
circ_seqs <- intersect(GenomeInfoDb::DEFAULT_CIRC_SEQS, seqlevels(gr) )
gr <- GenomicFeatures:::.tidy_seqinfo(gr, circ_seqs, NULL)
genome(seqinfo(gr))=genome
names(gr) <- gr$transcript_id
tximetaInfo <- list(version = packageVersion("tximeta"), importTime = Sys.time())
metadata <- list(tximetaInfo = tximetaInfo)
metadata$level <- "txp"
metadata$gtf <- GenomicFeatures:::.prepareGFFMetadata(annot, paste0(provider,"_v",version), organism, NA, NA, NULL)

assay.nms <- rownames(assay(txim,"counts"))
txps.missing <- !assay.nms %in% names(gr)
txps2 <- gr[rownames(assay(txim,"counts"))]
ucsc.genome <- tximeta:::genome2UCSC(genome)
try(seqinfo(txps2) <- Seqinfo(genome = ucsc.genome)[seqlevels(txps2)])

txomeInfo = list( index=gsub(annot,"",gtf.path),
                  source=provider,organism=organism,release=version,
                  genome =genome,
                  fasta= fasta,
                  gtf= annot,
                  sha256= paste0(provider,"_v",version,"_hash") )
rowRanges(txim) = txps2
metadata(txim)  = c(metadata(txim),metadata)
metadata(txim)$txomeInfo = txomeInfo

## add raw read counts from files
rcmat_files_list = list.files(".",pattern = "transcript_count_matrix",full.names = T)
rcmat_list = lapply(rcmat_files_list, read_csv,col_names = T)

for(x in rcmat_list) x = x[match(rownames(txim),x$transcript_id),]
rcmat = rcmat_list[[1]]
if(length(rcmat_list)>1){ 
    for(x in rcmat_list[-1]) rcmat = bind_cols(rcmat,y)[,colnames(txim)]
}
rcmatdf = as.data.frame(rcmat[,-1])
rownames(rcmatdf) = rcmat$transcript_id

assays(txim,withDimnames=F)$raw_counts = rcmatdf

names(assays(txim))[names(assays(txim))=="counts"] = "counts_length_normalized"
names(assays(txim))[names(assays(txim))=="abundance"] = "abundance_FPKM"

assays(txim)$abundance_TPM = sweep(assays(txim)$abundance_FPKM, 2,colSums(assays(txim)$abundance_FPKM),"/" )*10**6
assays(txim) <- SimpleList(counts=assay(txim,"raw_counts"),length=assay(txim,"length"),abundance_FPKM=assay(txim,"abundance_FPKM"),abundance_TPM=assay(txim,"abundance_TPM") )

eval(call("<-", as.name(paste0("transcript",suffix,".SE")),txim ))
save(list=paste0("transcript",suffix,".SE"),file = paste0("transcript",suffix,".SE.rda") )
rm(txim)
gc()

## 
# create Txome annotation for our annot
names(gr) <- gr$gene_id
assay.nms <- rownames(assay(txim.genes,"counts"))
genes.missing <- !assay.nms %in% names(gr)
genes2 <- gr[rownames(assay(txim.genes,"counts"))]
try(seqinfo(genes2) <- Seqinfo(genome = ucsc.genome)[seqlevels(genes2)])

metadata.genes = metadata
metadata.genes$level <- "gene"
rowRanges(txim.genes) = genes2
metadata(txim.genes)  = c(metadata(txim.genes),metadata.genes)
metadata(txim.genes)$txomeInfo = txomeInfo

## add raw read counts from files
names(assays(txim.genes))[names(assays(txim.genes))=="counts"] = "counts_length_normalized"

rcmat_files_list.g = list.files(".",pattern = "gene_count_matrix",full.names = T)
rcmat_list.g = lapply(rcmat_files_list.g, read_csv,col_names = T)

for(x in rcmat_list.g) x = x[match(rownames(txim.genes),x$gene_id),]
rcmat.g = rcmat_list.g[[1]]
if(length(rcmat_list.g)>1){ 
    for(x in rcmat_list.g[-1]) rcmat.g = bind_cols(rcmat.g,y)[,colnames(txim.genes)]
}
rcmatdf.g = as.data.frame(rcmat.g[,-1])
rownames(rcmatdf.g) = rcmat.g$gene_id

assays(txim.genes,withDimnames=F)$raw_counts = rcmatdf.g
names(assays(txim.genes))[names(assays(txim.genes))=="abundance"] = "abundance_FPKM"
assays(txim.genes)$abundance_TPM = sweep(assays(txim.genes)$abundance_FPKM, 2,colSums(assays(txim.genes)$abundance_FPKM),"/" )*10**6
assays(txim.genes) <- SimpleList(counts=assay(txim.genes,"raw_counts"),length=assay(txim.genes,"length"),abundance_FPKM=assay(txim.genes,"abundance_FPKM"),abundance_TPM=assay(txim.genes,"abundance_TPM") )

eval(call("<-", as.name(paste0("gene",suffix,".SE")),txim.genes ))
save(list=paste0("gene",suffix,".SE"),file = paste0("gene",suffix,".SE.rda") )
rm(txim.genes)
gc()