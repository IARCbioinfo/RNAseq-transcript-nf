require(ballgown)
samples = list.dirs(".")
samples = gsub("./", "", samples[samples!="."])
my.data = ballgown(samples=samples, dataDir=".", meas='all')
gdata = gexpr(my.data)
gdata = cbind(rownames(gdata),gdata)
colnames(gdata)[1] = "gene_id"
tdata = texpr(my.data)
tdata = cbind(transcriptNames(my.data),tdata)
colnames(tdata)[1] = "transcript_id"
write.csv(gdata,"gene_FPKM_matrix.csv",quote=F,row.names = F)
write.csv(tdata,"transcript_FPKM_matrix.csv",quote=F,row.names = F)
	
