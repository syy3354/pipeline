#enhancerPromoter.R

#plotting functions for enhancerPromoter code

#=========================================================
#===========================JOB PARAMETERS================
#=========================================================


args = commandArgs()
print(args)
geneTablePath = args[3]
outputFolder = args[4]
analysisName = args[5]
top = args[6]

#=========================================================
#===========================DEBUG SECTION=================
#=========================================================


# setwd('~/Dropbox/mycn_cyl')

# #needed arguments

# #gene table
# geneTablePath = './enhancerPromoter/HG19_P493-6_MYC_T24_-0_+0/HG19_P493-6_MYC_T24_-0_+0_GENE_TABLE.txt'

# #output folder
# outputFolder = './enhancerPromoter/HG19_P493-6_MYC_T24_-0_+0/output/'

# top = 5000


# analysisName = 'P493-6_MYC_T24'

#=========================================================
#========================HELPER FUNCTIONS=================
#=========================================================


plotContribution <- function(geneTable,analysisName,outputFolder,top=0,nBins=100){
	
	if(top == 0){
		top = nrow(geneTable)
		topString = 'all'
	}
	else if(top > nrow(geneTable)){
		top = nrow(geneTable)
		topString = 'all'			
	}else{
		topString = as.character(top)
	}

	#get the signal for tss, distal and total
	promoterSignal =geneTable[,2]
	enhancerSignal = geneTable[,3]
	totalSignal = promoterSignal + enhancerSignal

	#order by total signal
	totalOrder = order(totalSignal,decreasing=FALSE)[(length(totalSignal)-top+1):length(totalSignal)]


	#do a simpile moving average w/ a step size set by number of bins
	enhancerVector = c()
	promoterVector = c()
	
	i = 1
	stepSize = length(totalOrder)/nBins
	
	while(i < (length(totalOrder) - stepSize)){
	enhancerVector = c(enhancerVector,mean(enhancerSignal[totalOrder][i:(i+stepSize)]))
	promoterVector = c(promoterVector,mean(promoterSignal[totalOrder][i:(i+stepSize)]))
	i = i + stepSize/2
	}
	
	#check for crazy outliers
	if(max(totalSignal)/quantile(totalSignal,.99) > 2){
		yMax = as.numeric(quantile(totalSignal,.99)*1.5)
	}else{
		yMax = max(totalSignal)
	}
	
	#set the linewidth 
	linewidth = max(0.001,100/length(totalOrder))
	print(linewidth)
	
	plotPath = paste(outputFolder,analysisName,'_TOP_', topString,'_ORDERED_CONTRIB.pdf',sep='')
	pdf(file=plotPath,width = 11,height=8.5)
	par(mfrow=c(2,1))
	#we want to plot the bigger factor as the background
	if(max(enhancerSignal) > max(promoterVector)){
		plot(totalSignal[totalOrder],type='h',ylim =c(0,yMax),col='blue',xlab='All genes ranked by total signal',main=analysisName,lwd=0.25)
		lines(promoterSignal[totalOrder],type='h',col=rgb(1,0,0,0.3),lwd= linewidth)
	}else{
		plot(totalSignal[totalOrder],type='h',ylim =c(0,yMax),col='red',xlab='All genes ranked by total signal',main=analysisName,lwd=0.25)
		lines(enhancerSignal[totalOrder],type='h',col=rgb(0,0,1,0.3),lwd= linewidth)		
	}
	legend(0,0.75* yMax,c('promoter contribution','enhancer contribution'),fill=c('red','blue'))			

	enhancerContribVector = enhancerVector/(promoterVector+enhancerVector)*100
	enhancerContribVector[is.na(enhancerContribVector)] <- 0
	plot(1:length(enhancerContribVector),enhancerContribVector,type='p',col='blue',pch=16,ylab='% enhancer contribution to total signal',xaxt='n',xlab='',ylim =c(0,max(enhancerContribVector*1.15)))
	x=1:length(enhancerContribVector)
	lw1 = loess(enhancerContribVector ~x)
	lines(x,lw1$fitted,col='blue',lwd=2)
	dev.off()
	
	
}


runWaterfall <- function(geneTable,analysisName,outputFolder,top=0){
	if(top == 0){
		top = nrow(geneTable)
		topString = 'all'
	}
	else if(top > nrow(geneTable)){
		top = nrow(geneTable)
		topString = 'all'			
	}else{
		topString = as.character(top)
	}

	#get the signal for tss, distal and total
	promoterSignal =geneTable[,2]
	enhancerSignal = geneTable[,3]
	totalSignal = promoterSignal + enhancerSignal

	#order by total signal
	totalOrder = order(totalSignal,decreasing=FALSE)[(length(totalSignal)-top+1):length(totalSignal)]

	topTable = geneTable[totalOrder,]
	
	topEnhancerContrib = topTable[,3]/(topTable[,2]+topTable[,3])
	topEnhancerContrib_log = log2(topTable[,3]/topTable[,2])
	topEnhancerContrib_log[which(topEnhancerContrib_log== -Inf)] <- -8
	topEnhancerContrib_log[which(topEnhancerContrib_log== +Inf)] <- +8
	
	topEnhancerContrib_log[which(topEnhancerContrib_log < -8)] <- -8
	topEnhancerContrib_log[which(topEnhancerContrib_log > 8)] <- 8	
	
	enhancerContribOrder = order(topEnhancerContrib)
	topContribOrderedTable = cbind(topTable[enhancerContribOrder,],topEnhancerContrib[enhancerContribOrder], topEnhancerContrib_log[enhancerContribOrder])
	colnames(topContribOrderedTable)[4] = 'ENHANCER_CONTRIBUTION'
	colnames(topContribOrderedTable)[5] = 'ENHANCER_PROMOTER_RATIO'

	
	plotPath = paste(outputFolder,analysisName,'_TOP_', topString,'_WATERFALL.pdf',sep='')
	pdf(file=plotPath,width = 8,height =6)
	plot(1:top,1-topContribOrderedTable[,4],type='h',col='red',ylim=c(-1,1),ylab='Relative enhancer/promoter contribution',xlab=paste('Top',top,'genes as ranked by total signal',sep=' '))
	lines(-1*topContribOrderedTable[,4],type='h',col='blue')
	dev.off()

	#setting the color for the log2 waterfall
	colorSpectrum <- colorRampPalette(c("red","grey","grey","blue"))(100)

	#setting a color data range
	minValue <- -8
	maxValue <- 8
	color_cuts <- seq(minValue,maxValue,length=100)
	color_cuts <- c(min(topEnhancerContrib_log,na.rm=TRUE), color_cuts,max(topEnhancerContrib_log,na.rm=TRUE))


	#add one extra min color to even out sampling
	colorSpectrum <- c(colorSpectrum[1],colorSpectrum[1],colorSpectrum)

	colorVector = c()
	for(i in enhancerContribOrder){
		delta = topEnhancerContrib_log[i]
		color = colorSpectrum[max(which(color_cuts <= delta))]
		colorVector =c(colorVector,color)	
	}	
	plotPath = paste(outputFolder,analysisName,'_TOP_', topString,'_WATERFALL_LOG.pdf',sep='')
	pdf(file=plotPath,width = 8,height =6)
	plot(1:top,topContribOrderedTable[,5],type='h',col=colorVector,ylim=c(-8,8),ylab='log2 enhancer/promoter ratio',xlab=paste('Top',top,'genes as ranked by total signal',sep=' '))
	dev.off()	
	
	#writing out the rank ordered table
	tablePath = paste(outputFolder,analysisName,'_TOP_', topString,'_ORDERED.txt',sep='')
	write.table(topContribOrderedTable,file= tablePath,quote=FALSE,row.names=FALSE,sep='\t')
	
	#making the gct
	filename_gct= paste(outputFolder,analysisName,'_top_', topString,'.gct',sep='')
	gctMatrix =matrix(ncol=4,nrow=nrow(topContribOrderedTable))
	colnames(gctMatrix) = c('NAME','DESCRIPTION','PROMOTER','DISTAL')
	gctMatrix[,1]= as.character(topContribOrderedTable[,1])
	gctMatrix[,3]= topContribOrderedTable[,2]
	gctMatrix[,4]= topContribOrderedTable[,3]
	
	gctHeader = matrix(data='',ncol=4,nrow=3)
	gctHeader[1,1]='#1.2'
	gctHeader[2,1]=nrow(topContribOrderedTable)
	gctHeader[2,2]='2'
	gctHeader[3,]=c('NAME','DESCRIPTION','PROMOTER','DISTAL')
	gctCombined = rbind(gctHeader,gctMatrix)
	write.table(gctCombined,file=filename_gct,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)	
	
	
	#making the cls
	filename_cls= paste(outputFolder,analysisName,'_top_', topString,'.cls',sep='')
	clsTable = matrix(data='',ncol=3,nrow=3)
	clsTable[1,] =c(2,2,1)
	clsTable[2,1]=paste('#','PROMOTER',sep='')
	clsTable[2,2]='DISTAL'
	clsTable[3,1]='PROMOTER'
	clsTable[3,2]='DISTAL'
	write.table(clsTable,file=filename_cls,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)		
	
	
}
#========================================================
#========================DATA INPUT======================
#========================================================

geneTable = read.delim(geneTablePath,sep='\t')




#========================================================
#=================PLOTTING CONTRIBUTION==================
#========================================================

#for all
plotContribution(geneTable,analysisName,outputFolder)
runWaterfall(geneTable,analysisName,outputFolder)


#top 5000
plotContribution(geneTable,analysisName,outputFolder,5000) 
runWaterfall(geneTable,analysisName,outputFolder,5000)

