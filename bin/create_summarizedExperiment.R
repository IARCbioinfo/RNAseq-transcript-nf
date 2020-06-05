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
library(GenomicFeatures)

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
readlengths_rcmat_list = sapply( rcmat_files_list, function(x) as.numeric(strsplit(strsplit(rev( strsplit(x,"_")[[1]] )[1],"l")[[1]][2],".csv")[[1]][1] ) )
samples_rcmat_list = lapply(rcmat_list, function(x) colnames(x)[-1])
readlengths_txim = as.data.frame(matrix(rep(NA,ncol(txim)),nrow=1 ))
colnames(readlengths_txim) = colnames(txim)
for(i in 1:length(rcmat_list)) readlengths_txim[1,samples_rcmat_list[[i]]] = readlengths_rcmat_list[i]

rcmat = rcmat_list[[1]]
if(length(rcmat_list)>1){ 
    for(i in 2:length(rcmat_list)) rcmat = left_join(rcmat,rcmat_list[[i]],by="transcript_id")
}
rcmat = rcmat[match(rownames(txim),rcmat$transcript_id),c("transcript_id",colnames(txim))]

write_csv(rcmat,path=paste0("transcript_count_matrix",suffix,".csv"))
rcmatdf = as.data.frame(rcmat[,-1])
rownames(rcmatdf) = rcmat$transcript_id

colData(txim)$readlength = readlengths_txim[1,]
assays(txim,withDimnames=F)$raw_counts = rcmatdf

names(assays(txim))[names(assays(txim))=="counts"] = "counts_length_normalized"
names(assays(txim))[names(assays(txim))=="abundance"] = "abundance_FPKM"

assays(txim)$abundance_TPM = sweep(assays(txim)$abundance_FPKM, 2,colSums(assays(txim)$abundance_FPKM),"/" )*10**6
assays(txim) <- SimpleList(counts=assay(txim,"raw_counts"),length=assay(txim,"length"),abundance_FPKM=assay(txim,"abundance_FPKM"),abundance_TPM=assay(txim,"abundance_TPM") )

# add exon information
#exons <- exonsBy(mytxdb, by="tx",use.names=T) #tximeta:::getRanges(txdb = mytxdb, txomeInfo = txomeInfo, type = "exon")
exons <- split(gr[gr$type=="exon",c("exon_id" ,"exon_name", "exon_number")[c("exon_id" ,"exon_name", "exon_number")%in%colnames(elementMetadata(gr))]],
                names(gr)[gr$type=="exon"])
exons <- exons[names(txim)]
if (all(is.na(seqlengths(exons)))) {
    seqinfo(exons) <- seqinfo(txim)
}
exons <- exons[rownames(txim)]
stopifnot(all(rownames(txim) == names(exons)))
mcols(exons) <- mcols(txim)
rowRanges(txim) <- exons

eval(call("<-", as.name(paste0("transcript",suffix,".SE")),txim ))
save(list=paste0("transcript",suffix,".SE"),file = paste0("transcript",suffix,".SE.rda") )
rm(txim)
gc()

## 
gr.g = gr[gr$type=="exon",]
names(gr.g) <- gr.g$gene_id
exons.g <- split(gr.g,names(gr.g) ) #groupGRangesBy(gr.g[gr.g$type=="exon",])
exons.g <- exons.g[names(txim.genes)]

genes2 = elementMetadata(unlist(exons.g))
rownames(genes2) = genes2$gene_id
genes2 = genes2[cumsum(elementNROWS(exons.g)),]
genes2$type = "gene"

assay.nms <- rownames(assay(txim.genes,"counts"))
genes.missing <- !assay.nms %in% names(gr.g)
genes2 <- genes2[rownames(assay(txim.genes,"counts")),]
try(seqinfo(genes2) <- Seqinfo(genome = ucsc.genome)[seqlevels(genes2)])

metadata.genes = metadata
metadata.genes$level <- "gene"
rowData(txim.genes) = genes2
metadata(txim.genes)  = c(metadata(txim.genes),metadata.genes)
metadata(txim.genes)$txomeInfo = txomeInfo

## add raw read counts from files
names(assays(txim.genes))[names(assays(txim.genes))=="counts"] = "counts_length_normalized"

rcmat_files_list.g = list.files(".",pattern = "gene_count_matrix",full.names = T)
rcmat_list.g = lapply(rcmat_files_list.g, read_csv,col_names = T)

rcmat.g = rcmat_list.g[[1]]
if(length(rcmat_list.g)>1){ 
    for(i in 2:length(rcmat_list.g)) rcmat.g = left_join(rcmat.g,rcmat_list.g[[i]],by="gene_id")
}
rcmat.g = rcmat.g[match(rownames(txim.genes),rcmat.g$gene_id),c("gene_id",colnames(txim.genes))]
write_csv(rcmat.g,path=paste0("gene_count_matrix",suffix,".csv"))
rcmatdf.g = as.data.frame(rcmat.g[,-1])
rownames(rcmatdf.g) = rcmat.g$gene_id

colData(txim.genes)$readlength = readlengths_txim[1,]

assays(txim.genes,withDimnames=F)$raw_counts = rcmatdf.g
names(assays(txim.genes))[names(assays(txim.genes))=="abundance"] = "abundance_FPKM"
assays(txim.genes)$abundance_TPM = sweep(assays(txim.genes)$abundance_FPKM, 2,colSums(assays(txim.genes)$abundance_FPKM),"/" )*10**6

### check if still SimpleList in recent SummarizedExperiment versions
assays(txim.genes) <- SimpleList(counts=assay(txim.genes,"raw_counts"),length=assay(txim.genes,"length"),abundance_FPKM=assay(txim.genes,"abundance_FPKM"),abundance_TPM=assay(txim.genes,"abundance_TPM") )

if (all(is.na(seqlengths(exons.g)))) {
    seqinfo(exons.g) <- seqinfo(exons)
}
exons.g <- exons.g[rownames(txim.genes)]
stopifnot(all(rownames(txim.genes) == names(exons.g)))
mcols(exons.g) <- mcols(txim.genes)
rowRanges(txim.genes) <- exons.g

eval(call("<-", as.name(paste0("gene",suffix,".SE")),txim.genes ))
save(list=paste0("gene",suffix,".SE"),file = paste0("gene",suffix,".SE.rda") )
rm(txim.genes)
gc()