#enhancerPromoter.R

#plotting functions for enhancerPromoter code

#=========================================================
#===========================JOB PARAMETERS================
#=========================================================


args = commandArgs()
print(args)
promoterTablePath = args[3]
distalTablePath = args[4]
outputFolder = args[5]
analysisName = args[6]
top = args[7]

#=========================================================
#===========================DEBUG SECTION=================
#=========================================================



#distalTablePath = './gsea/NB_MYCN_CONSERVED_FILTERED_top_2000.Gsea.1459185088639/gsea_report_for_DISTAL_1459185088639.xls'
#promoterTablePath = './gsea/NB_MYCN_CONSERVED_FILTERED_top_2000.Gsea.1459185088639/gsea_report_for_PROMOTER_1459185088639.xls'

#outputFolder = './'
#analysisName = 'NB_MYCN_CONSERVED_FILTERED'
#top = '2000'
#========================================================
#========================DATA INPUT======================
#========================================================

distalTable = read.delim(distalTablePath)
promoterTable = read.delim(promoterTablePath)



#========================================================
#====================MAKING NES TABLE====================
#========================================================

nes_table = rbind(promoterTable[,c(1,4,6,8,10)],distalTable[,c(1,4,6,8,10)])


nesOrder = order(nes_table$NES)
nes_table = nes_table[nesOrder,]

#now write to disk

nes_table_path = paste(outputFolder,analysisName,'_top_',top,'_nes.txt',sep='')

write.table(nes_table,nes_table_path,quote=FALSE,row.names=FALSE,sep='\t')


#========================================================
#======================PLOTTING NES======================
#========================================================


nes_plot_path = paste(outputFolder,analysisName,'_top_',top,'_nes.pdf',sep='')
pdf(file= nes_plot_path,width = 5,height =5)
plot(nes_table[,3],nes_table[,4],pch=16,cex=1,col='grey',ylab='FDR',xlab='NES')
abline(h=0.1,lty=2)

dev.off()