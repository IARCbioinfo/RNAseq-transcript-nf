require(ballgown)
my.data = ballgown(dataDir=".", meas='all')
gdata = gexpr(my.data)
tdata = texpr(my.data)
rownames(tdata) = transcriptNames(my.data)
write.csv(gdata,"gene_FPKM_matrix.csv",quote=F)
write.csv(tdata,"transcript_FPKM_matrix.csv",quote=F)
	