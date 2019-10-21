args = commandArgs(trailingOnly=TRUE)
suffix=args[1]
if(is.na(suffix)) suffix = ""

require(ballgown)

#load samples
samples = list.dirs(".")
samples = gsub("./", "", samples[samples!="."])
bg = ballgown(samples=samples, dataDir=".", meas='all')
gdata = gexpr(bg)
gdata = cbind(rownames(gdata),gdata)
colnames(gdata)[1] = "gene_id"
tdata = texpr(bg)
tdata = cbind(transcriptNames(bg),tdata)
colnames(tdata)[1] = "transcript_id"

#write tables and R object
write.csv(gdata,paste0("gene_FPKM_matrix",suffix,".csv"),quote=F,row.names = F)
write.csv(tdata,paste0("transcript_FPKM_matrix",suffix,".csv"),quote=F,row.names = F)
save(bg,file=paste0("bg",suffix,".rda") )
