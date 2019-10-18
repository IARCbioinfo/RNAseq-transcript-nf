require(ballgown)
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
write.csv(gdata,"gene_FPKM_matrix.csv",quote=F,row.names = F)
write.csv(tdata,"transcript_FPKM_matrix.csv",quote=F,row.names = F)
save(bg,file="bg.rda")
