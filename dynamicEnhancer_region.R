#dynamicEnhancer_region_plot.R

#makes line and bar plots for all regions in a dynamic enhancer meta analysis 

# The MIT License (MIT)

# Copyright (c) 2016 Charles Lin

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

#enhancerFile ='mergeTest/EC_BRD4_CON_ROSE/HG18_EC_MERGED_SUPERS_-0_+0_0KB_STITCHED_ENHANCER_DELTA_ENHANCER_TO_GENE_100KB_RANK.txt'

#nSuper1 = 347
#nSuper2 = 271
#===========================================================
#===============READING IN ARGUMNETS========================
#===========================================================

args <- commandArgs()

print(args)

normFile = args[6]


name1 = args[7]
name2 = args[8]

group1_length = as.numeric(args[9])
group2_length = as.numeric(args[10])


#===========================================================
#=======================TESTING=============================
#===========================================================

# normFile = '/Volumes/grail/projects/gbm/bordo/gbm/dynamic_meta_rose/GBM_H2S-TPC_H3K27AC_ROSE/HG19_GROUP1_GROUP2_merged_MERGED_REGIONS_-0_+0_0KB_STITCHED_ENHANCER_REGION_MAP_NORM.txt'

# name1 = 'GROUP1'
# name2 = 'GROUP2'

# group1_length = 4
# group2_length = 5


#===========================================================
#===================HELPER FUNCTIONS========================
#===========================================================

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
        if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
        stop("vectors must be same length")
        arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
        }




rankPlotRegion <- function(diffTable,rankMatrixInverted_1, rankMatrixInverted_2,name1,name2,maxRank=100){
	
	#uses 1/log(rank) scaling for enhancers where rank is from 1-nEnancers with 1 being best enhancer
	#transforms a matrix where rank is given in increasing order (e.g. enhancer rank 1000 is higher than 500)

	group1_length = ncol(rankMatrixInverted_1)
	group2_length = ncol(rankMatrixInverted_2)
	#set a maxLog y limit
	scaleBaseStart =7
	yLimit = nrow(rankMatrixInverted_1)
	cut_off = 1000
	

	#set a yTop and yBottom
	yTop = 1/log(maxRank)
	yBottom = 1/log(yLimit)	
	
	colorVector = c(rep('grey', group1_length),rep('red', group2_length))	

	for(i in 1:nrow(diffTable)){

		region_id = as.character(diffTable[i,1])
		location_string = paste(diffTable[i,2],':', diffTable[i,3],' ', diffTable[i,4],sep='')
		region_row = which(rownames(rankMatrixInverted_1) == region_id)		
		rankVector = c(rankMatrixInverted_1[region_row,],rankMatrixInverted_2[region_row,])
	
		par(mai=c(3,1,1,0.25))
		layout(matrix(c(1,1,2),nrow=1))
		plot(0,0,xlim = c(0.25, group1_length + group2_length +.5),ylim=c(yBottom,yTop),cex=0,yaxt='n',xaxt='n',xlab='',ylab='Region rank in sample',main = c(region_id,location_string))
		
		#y axis
		yTicks = c(maxRank,100,500,1000)
		yTicks = 1/log(yTicks)
		axis(2,yTicks, c(paste('<',as.character(maxRank)),100,500,1000),las=2)
		
		xTicks = c(.5, group1_length +0.5, group1_length + group2_length +0.5)
		xLabTicks = c(mean(xTicks[1:2]),mean(xTicks[2:3]))
		axis(1,xTicks,labels=FALSE)
		#axis(1,xLabTicks,c(name1,name2),tick=FALSE)
		
		axis(1,1:length(colorVector),c(colnames(rankMatrixInverted_1),colnames(rankMatrixInverted_2)),las=3,tick=FALSE,hadj=1,col=colorVector)

		alphaScale = 50
		polygon(c(.5,.5,xTicks[2],xTicks[2]),c(yBottom,yTop,yTop,yBottom),col=rgb(200,200,200,alphaScale,maxColorValue=255),border=NA)	
		polygon(c(xTicks[2],xTicks[2],xTicks[3],xTicks[3]),c(yBottom,yTop,yTop,yBottom),col=rgb(200,0,0,alphaScale,maxColorValue=255),border=NA)	
		
		lines(1:length(rankVector),1/log(rankVector),col='black',type= 'l',cex=0)
		abline(h=1/log(cut_off),col='grey',lty=2)
		topPoints = which(rankVector <= cut_off)	
		points(topPoints,1/log(rankVector[topPoints]),col=colorVector[topPoints],cex=1.5,pch=19)	
		
		
		#now do the dot plot
		#get the values for group1
		group1_sig = as.numeric(diffTable[i,7:(7+group1_length-1)])
		group2_sig = as.numeric(diffTable[i,(7+group1_length):(ncol(diffTable)-5)])
		
		xRange = c(0,2)
		yRange = c(0,1.25*max(group1_sig,group2_sig))
		plot(0,0,ylim = yRange,xlim = xRange,cex=0,xaxt='n',ylab = 'median fold region signal',xlab='',las=2)
		axis(1,c(0.5,1.5),labels = c(name1,name2),las=2)
		xVector_1 = jitter(rep(0.5,length(group1_sig)),amount =0.1)
		xVector_2 = jitter(rep(1.5,length(group2_sig)),amount =0.1)

		segments(.4,mean(group1_sig),.6,mean(group1_sig),lwd=4,col = 'black')
		segments(1.4,mean(group2_sig),1.6,mean(group2_sig),lwd=4,col = 'red')

		points(xVector_1,group1_sig,col = rgb(0.5,0.5,0.5,0.5),pch=16,cex=1.5)
		points(xVector_2,group2_sig,col = rgb(1,0,0,0.5),pch=16,cex=1.5)
		

	}				
	
	

}
	
	
    
#===========================================================
#==================GETTING IN TABLES========================
#===========================================================

normTable = read.delim(normFile)


#first make a normMatrix for each group

normMatrix_1 = as.matrix(normTable[,7:(7+group1_length-1)])
normMatrix_2 = as.matrix(normTable[,(7+group1_length):ncol(normTable)])

rownames(normMatrix_1) = as.character(normTable$REGION_ID)
rownames(normMatrix_2) = as.character(normTable$REGION_ID)

#make rank matrices 
#higher rank = higher signal
rankMatrix_1 = apply(normMatrix_1,2,rank)
rankMatrix_2 = apply(normMatrix_2,2,rank)

maxRank = 100
rankMatrixInverted_1 = nrow(rankMatrix_1)-rankMatrix_1+1
rankMatrixInverted_1[which(rankMatrixInverted_1 <maxRank)] <- maxRank  #sets a ceiling on rank

rankMatrixInverted_2 = nrow(rankMatrix_2)-rankMatrix_2+1
rankMatrixInverted_2[which(rankMatrixInverted_2 <maxRank)] <- maxRank  #sets a ceiling on rank


#===========================================================
#==============FINDING DIFFERENTIAL REGIONS=================
#===========================================================
#now let's do some actual differential analysis
#need guys w/ fold change > log2(1) AND fdr less than 0.05

p_vector_wilcox = c()

fold_vector = c()
for(i in 1:nrow(normMatrix_1)){
	
	#t-test
	p_wilcox = wilcox.test(normMatrix_1[i,],normMatrix_2[i,])$p.value
	p_vector_wilcox = c(p_vector_wilcox,p_wilcox)
	
	
	fold = log2(mean(normMatrix_2[i,])/mean(normMatrix_1[i,]))
	fold_vector = c(fold_vector,fold)
}

#sig_rows = which(p.adjust(p_vector_wilcox,'fdr')<= 0.05)
sig_rows = which(p_vector_wilcox <= 0.05)

fold_rows = which(abs(fold_vector) >= 1 )
diff_rows = intersect(sig_rows,fold_rows)

fold_rank = length(fold_vector) - rank(fold_vector)
is_diff = rep(0,length(fold_vector))
is_diff[diff_rows] <- 1

newTable = cbind(normTable,fold_vector,p_vector_wilcox,p.adjust(p_vector_wilcox,'fdr'),fold_rank,is_diff)
colnames(newTable)[(ncol(normTable)+1):ncol(newTable)] = c(paste(name2,'_vs_',name1,'_LOG2_FOLD',sep=''),'P_VALUE_WILCOX','FDR_VALUE_WILCOX','RANK','DIFF_0.05')

#now make a diff table that is ranked by fold change

diffTable = newTable[diff_rows,]
diffOrder = order(fold_vector[diff_rows],decreasing=TRUE)
diffTable = diffTable[diffOrder,]

#round the tables to make them less sad


newTable[,7:ncol(newTable)] = round(newTable[,7:ncol(newTable)],4)
diffTable[,7:ncol(diffTable)] = round(diffTable[,7:ncol(diffTable)],4)


#now we need to write these two tables out to disk
region_out = gsub('REGION_MAP_NORM.txt','REGION_STATS.txt',normFile)
diff_out = gsub('REGION_MAP_NORM.txt','REGION_STATS_DIFF.txt',normFile)

write.table(newTable,file=region_out,quote=FALSE,row.names=FALSE,col.names=TRUE,sep='\t')
write.table(diffTable,file= diff_out,quote=FALSE,row.names=FALSE,col.names=TRUE,sep='\t')


#===========================================================
#============RANK PLOTS OF DIFFERENTIAL REGIONS=============
#===========================================================


#first on gained regions


plot_width = .7*(ncol(rankMatrixInverted_1) + ncol(rankMatrixInverted_2))+2

#SETTING THE OUTPUTS
gainedRows = which(diffTable[,(ncol(diffTable)-4)] >= 1)
rank_gained_out = gsub('REGION_MAP_NORM.txt','REGION_GAINED.pdf',normFile)

lostRows = which(diffTable[,(ncol(diffTable)-4)] <= -1)
rank_lost_out = gsub('REGION_MAP_NORM.txt','REGION_LOST.pdf',normFile)


unchangedRows = which(newTable[,ncol(newTable)] ==0)
unchanged_out = gsub('REGION_MAP_NORM.txt','REGION_UNCHANGED.pdf',normFile)




pdf(file = rank_gained_out,width = plot_width ,height =8.5)
rankPlotRegion(diffTable[gainedRows,],rankMatrixInverted_1, rankMatrixInverted_2,name1,name2,maxRank=100)
dev.off()

pdf(file = rank_lost_out,width = plot_width ,height =8.5)
rankPlotRegion(diffTable[lostRows,],rankMatrixInverted_1, rankMatrixInverted_2,name1,name2,maxRank=100)
dev.off()

pdf(file = unchanged_out,width = plot_width ,height =8.5)
rankPlotRegion(newTable[unchangedRows,],rankMatrixInverted_1, rankMatrixInverted_2,name1,name2,maxRank=100)
dev.off()



