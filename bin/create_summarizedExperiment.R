args = commandArgs(trailingOnly=TRUE)
annot=args[1]
provider=args[2]
version=args[3]
genome=args[4]
organism=args[5]
gtf.path=args[6]

samples = list.dirs(".")
samples = gsub("./", "", samples[samples!="."])

# transcriptome R file creation
library(tximeta)
library(tximport)
library(readr)
library(dplyr)
#require(SummarizedExperiment)

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
gr <- GenomicFeatures:::.tidy_seqinfo(gr, GenomicFeatures:::DEFAULT_CIRC_SEQS, NULL)
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
                  source=provider,organism=organism,release=31,
                  genome =genome,
                  fasta= "",#paste0(gsub(annot,"",gtf.path),"ref_genome.fa"),
                  gtf= annot,
                  sha256= paste0(provider,"_v",version,"_hash") )
rowRanges(txim) = txps2
metadata(txim)  = c(metadata(txim),metadata)
metadata(txim)$txomeInfo = txomeInfo

## add raw read counts from files
rcmat_files_list = list.files(".",pattern = "transcript_count_matrix",full.names = T)
rcmat_list = lapply(rcmat_files_list, read_csv,col_names = T)
#rcmat_l75 = read_csv("transcript_count_matrix_l75.csv",col_names = T)

for(x in rcmat_list) x = x[match(rownames(txim),x$transcript_id),]
#rcmat_l75  = rcmat_l75[match(rownames(txim),rcmat_l75$transcript_id),]
rcmat = rcmat_list[[1]]
if(length(rcmat)>1){ 
    for(x in rcmat_list[-1]) rcmat = bind_cols(rcmat,y)[,colnames(txim)]
}

assays(txim)$raw_counts = rcmat

names(assays(txim))[names(assays(txim))=="counts"] = "counts_length_normalized"
names(assays(txim))[names(assays(txim))=="abundance"] = "abundance_FPKM"

assays(txim)$abundance_TPM = sweep(assays(txim)$abundance_FPKM, 2,colSums(assays(txim)$abundance_FPKM),"/" )*10**6
assays(txim) <- SimpleList(counts=assay(txim,"counts"),length=assay(txim,"length"),abundance_FPKM=assay(txim,"abundance_FPKM"),abundance_TPM=assay(txim,"abundance_TPM") )

transcript.SE = txim
rm(txim)
gc()

save(transcript.SE,file = "transcript.rda")

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

for(x in rcmat_list.g) x = x[match(rownames(txim),x$gene_id),]
#rcmat_l75  = rcmat_l75[match(rownames(txim),rcmat_l75$transcript_id),]
rcmat.g = rcmat_list.g[[1]]
if(length(rcmat.g)>1){ 
    for(x in rcmat_list.g[-1]) rcmat.g = bind_cols(rcmat.g,y)[,colnames(txim.genes)]
}


#rcmat_l100.g = read_csv("gene_count_matrix_l100.csv",col_names = T)
#rcmat_l100.g$gene_id = sapply(rcmat_l100.g$gene_id, function(x) strsplit(x,"\\|")[[1]][1] )
#rcmat_l75.g  = read_csv("gene_count_matrix_l75.csv",col_names = T)
#rcmat_l75.g$gene_id = sapply(rcmat_l75.g$gene_id, function(x) strsplit(x,"\\|")[[1]][1] )
#rcmat_l100.g = rcmat_l100.g[match(rownames(txim.genes),rcmat_l100.g$gene_id),]
#rcmat_l75.g  = rcmat_l75.g[match(rownames(txim.genes),rcmat_l75.g$gene_id),]
#rcmat.g = bind_cols(rcmat_l100.g,rcmat_l75.g)[,colnames(txim.genes)]

assays(txim.genes)$raw_counts = rcmat.g
names(assays(txim.genes))[names(assays(txim.genes))=="abundance"] = "abundance_FPKM"
assays(txim.genes)$abundance_TPM = sweep(assays(txim.genes)$abundance_FPKM, 2,colSums(assays(txim.genes)$abundance_FPKM),"/" )*10**6
assays(txim.genes) <- SimpleList(counts=assay(txim.genes,"counts"),length=assay(txim,"length"),abundance_FPKM=assay(txim.genes,"abundance_FPKM"),abundance_TPM=assay(txim.genes,"abundance_TPM") )

gene.SE = txim.genes
rm(txim.genes)
gc()

save(gene.SE,file = "gene.SE.rda")

