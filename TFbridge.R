##Check whether the required packages are installed or not
m.pkgs <- installed.packages()
pkgs <- rownames(m.pkgs)

required.pkgs <- c("optparse", "preprocessCore", "MASS", "MXM", "parallel")
if(!(required.pkgs[1] %in% pkgs)){
    stop(sprintf("Please install %s CRAN package before proceeding", 
    	required.pkgs[1]))
}
if(!(required.pkgs[2] %in% pkgs)){
    stop(sprintf("Please install %s Bioconductor package before proceeding", 
    	required.pkgs[2]))
}
if(!(required.pkgs[3] %in% pkgs)){
    stop(sprintf("Please install %s CRAN package before proceeding", required.pkgs[3]))
}
if(!(required.pkgs[4] %in% pkgs)){
    stop(sprintf("Please install %s CRAN package before proceeding", required.pkgs[4]))
}
if(!(required.pkgs[5] %in% pkgs)){
    stop(sprintf("Please install %s CRAN package before proceeding", required.pkgs[5]))
}

##Quantile normalize the features
normalizeMatrix <- function(m, doLog2 = T){

  if(doLog2){
    ##Log2 of the values
    cu <- log2(m + 1)
  } else {
    cu <- m
  }

  ##Quantile normalize the values
  cu.qq <- normalize.quantiles(cu)
  colnames(cu.qq) <- colnames(m)
  rownames(cu.qq) <- rownames(m)

  return(cu.qq)

}

##This function returns the names of the genes that are expressed in atleast
##given proportion of samples by a given cutoff
getExpressedGenes <- function(df.rpkm, samples = NULL, present.prop = 0.5, rpkm.cutoff = 1,
  need.symbol = T, need.values = F){

  m.rpkm.grp <- df.rpkm[ ,samples]
  rownames(m.rpkm.grp) <- df.rpkm$GENE_SYMBOL
  exprs.logical <- apply(m.rpkm.grp, 1, function(x) sum(x > 1) > (present.prop * length(x)))
  m.rpkm.grp.exprs <- m.rpkm.grp[exprs.logical, ]

  if(need.values){
    return(m.rpkm.grp.exprs)
    } else{
      return(rownames(m.rpkm.grp.exprs))
    }

}

##This function gives all the predicted targets from CRC for all the given samples
getALLTargetGenes <- function(df.sample.info = NULL, samples = NULL){
  all.targets <- c()
  for(sample in samples){
    print(sample)
    ##location of the edge table
    edge.table.file.location <- as.character(subset(df.sample.info, SAMPLE_NAME == sample)$FILE)
    edge.table.file <- file.path(edge.table.file.location, sample, sprintf("%s_EDGE_TABLE.txt", sample))
        if(!file.exists(edge.table.file)){
          stop(sprintf("%s file does not exist", edge.table.file))
        }   
    df.edg <- read.delim(edge.table.file, sep = "\t", header = T, stringsAsFactors = F)
    sample.targets <- unique(df.edg$TARGET)
    all.targets <- unique(c(all.targets, sample.targets))
  }
  return(all.targets)
}

##This function prepares a matrix of all the TFs and its targets
##It works as union of all samples
##Column - TFs
##Rows - Target genes
getMotifHitMatrix <- function(df.sample.info = NULL, samples = NULL){

  df.all.edges <- data.frame()
  for(sample in samples){
    ##location of the edge table
    edge.table.file.location <- as.character(subset(df.sample.info, SAMPLE_NAME == sample)$FILE)
    edge.table.file <- file.path(edge.table.file.location, sample, sprintf("%s_EDGE_TABLE.txt", sample))
        if(!file.exists(edge.table.file)){
          stop(sprintf("%s file does not exist", edge.table.file))
        }
    df.edg <- read.delim(edge.table.file, sep = "\t", header = T, stringsAsFactors = F)
    df.all.edges <- rbind(df.all.edges, df.edg)
  }

  all.targets <- sort(unique(df.all.edges$TARGET))
  all.sources <- sort(unique(df.all.edges$SOURCE))

  m.motif.hit <- as.matrix(do.call(cbind, sapply(all.sources, function(x){
    this.target <- subset(df.all.edges, SOURCE == x)$TARGET
    vec0 <- mat.or.vec(nr = 1, nc = length(all.targets))
    names(vec0) <- all.targets
    vec0[this.target] <- 1
    return(vec0)
    }, simplify = F)))

  return(m.motif.hit)
}

##This will fit the Ridge for every target gene
fitRidge <- function(target.gene, df.rpkm, samples, m.motif.hit, Nrand = 10000){

  ##Get the exprssed gene matrix and quantile normalize it
  m.rpkm.grp.exprs <- getExpressedGenes(df.rpkm, samples = samples, need.values = T)
  rownames(m.rpkm.grp.exprs) <- toupper(rownames(m.rpkm.grp.exprs))
  m.rpkm.grp.exprs <- normalizeMatrix(as.matrix(m.rpkm.grp.exprs))

  ##Find the TFs that bind into the target gene
  vec.motif.hit <- m.motif.hit[target.gene, ]
  binding.tfs <- names(which(vec.motif.hit == 1))
  binding.tfs <- setdiff(binding.tfs, target.gene)
  binding.tfs <- intersect(rownames(m.rpkm.grp.exprs), binding.tfs)

  ##Get the expression value of all the binding Tfs and 
  ##the corresponding target gene expression for every sample
  m.X <- {}
  vec.Y <- c()
  for(every.sample in samples){
    exprs.vals <- m.rpkm.grp.exprs[binding.tfs, every.sample]
    m.X <- rbind(m.X, exprs.vals)
    vec.Y <- c(vec.Y, m.rpkm.grp.exprs[target.gene, every.sample])
  }
  m.data <- cbind(m.X, y = vec.Y)
  m.data <- scale(m.data, center = T, scale = F)
  rownames(m.data) <- NULL
  df.data <- as.data.frame(m.data)

  ##Learn the lambda by 10 fold cross validation
  lambda.val <- ridgereg.cv(m.data[ ,ncol(m.data)], m.data[ ,-ncol(m.data)], K = 10, lambda = seq(0, 100, by = 0.01), 
   auto = TRUE, seed = FALSE, ncores = 10, mat = NULL)$lambda
  ##Use the lambda for fitting the ridge
  model <- lm.ridge(y ~ ., df.data, lambda = lambda.val)
  pvals <- permuteAndGetCoefficients(df.data, N = Nrand, model$coef, lambda.val)
  coef <- model$coef
  return(list(pvalue = pvals, coefficent = coef))

}

##Permute the Y and refit the ridge
##This will give empirical distribution of the coefficients
##Assess if the coefficents learnt from the permuted model are greater than
##real coef
##Redo the above methodolgy N number of times and find the proportion of times
##the random coefficents are greater than real coefficent
permuteAndGetCoefficients <- function(df.data, N = 10000, real.coef, lambda.val){

  m.all.permuted.coeff <- matrix()
  m.coeff <- do.call(rbind, mclapply(1:N, function(x){
    df.permuted.data <- df.data
    df.permuted.data$y <- sample(df.data$y, size = length(df.data$y))
    model <- lm.ridge(y ~ ., df.permuted.data, lambda = lambda.val)
    model$coef
  }, mc.cores = 20))
  m.coeff <- abs(m.coeff)
  all.pvalue <- c()
  for(idx in 1:ncol(m.coeff)){
    pvalue <- length(which(m.coeff[ ,idx] > abs(real.coef[idx])))/N
    all.pvalue <- c(all.pvalue, pvalue)
  }
  names(all.pvalue) <- setdiff(colnames(df.data), "y")
  return(all.pvalue)
}


##For every target gene find the TFs are important
##A TF that does not bind gets a value 0
##A TF with pval <= pval.cutoff get a value 1
##A TF with pval > pval.cutoff get a value 0
##At the end it returns a matrix, target genes in rows and 
#TFs in columns, with the important TFs with value 1
##All the target genes that had no significant TF or all the 
getPVALMatrix <- function(ls.motif.pvals, m.motif.hit, pval.cutoff = 0.05){
  
  target.genes <- names(ls.motif.pvals)
  mat <- do.call(rbind, sapply(target.genes, function(x){
    motif.pval <- ls.motif.pvals[[x]]$pvalue
    vec.motif.hit <- m.motif.hit[x, ]
    vec.motif.hit[vec.motif.hit == 0] <- 1
    motif.pval[motif.pval > pval.cutoff] <- 1
    motif.pval[motif.pval <= pval.cutoff] <- 0.1
    vec.motif.hit[names(motif.pval)] <- motif.pval
    vec.motif.hit1 <- -log10(vec.motif.hit)
    }, simplify = F))
  #rm.all0.cols <- apply(mat, 2, function(x) !all(x == 0))
  #rm.all0.rows <- apply(mat, 1, function(x) !all(x == 0))
  #mat <- mat[rm.all0.rows ,rm.all0.cols]

  return(mat)

}

##This function returns the estimated coefficient and the significant
##of the interaction beween the TFs and the target genes
getAllInteractions <- function(ls.motif.pvals){

	target.genes <- names(ls.motif.pvals)
	df.interactions <- do.call(rbind, sapply(target.genes, function(x){
	    pvalue <- ls.motif.pvals[[x]]$pvalue
	    coefficient <- ls.motif.pvals[[x]]$coefficent
	    all.tfs <- unique(names(pvalue), names(coefficient))
	    pvalue.ordered <- pvalue[all.tfs]
	    coefficient.ordered <- coefficient[all.tfs]
	    target.gene.name <- rep(x, length(all.tfs))
	    df <- data.frame(tf = all.tfs, target_gene = target.gene.name, pvalue = pvalue.ordered, coefficient = coefficient.ordered)
	    }, simplify = F))
	rownames(df.interactions) <- NULL
	return(df.interactions)
	
}

########################################################################################################
##Load the required packages
suppressMessages(library(optparse))
suppressMessages(library(preprocessCore))
suppressMessages(library(parallel))
suppressMessages(library(MASS))
suppressMessages(library(MXM))

##List of the command line arguments
option_list = list(
  make_option(c("-f", "--gfile"), type="character", default=NULL, 
              help="Gene expression file name (RPKM/FPKM). One of the columns in the file should GENE_SYMBOL.
                The gene symbols should be consistent with ones used in the edge table. The other columns will have sample
                names as headers."),
	make_option(c("-d", "--datatable"), type="character", default=NULL, 
              help="Data table file name with sample information. This should be a tab delimited 
                file with atleast two columns. Column1:header SAMPLE_NAME which gives the name of the samples.
                These names should be consistent with the sample names used in the gene expression file.
                Column2: header FILE that gives the location of the directory where
                the edge table file is located (basically looking for your CRC output)"),
	make_option(c("-r", "--Nrandom"), type="integer", default=NULL, 
              help="Number of randomizations for fitting ridge regression for obtaining 
                empirical distribution over the TF coefficients. If no value is given then this number will be
                calculated based on the number of samples in the dataset. If the number of possible randomizations
                is more than 10,000 then it does 10,000 maximum randomizations"),
	make_option(c("-e", "--expressioncutoff"), type="double", default=1, 
              help="Gene expression cutoff. 
                Default is 1 (RPKM/FPKM)"),
	make_option(c("-s", "--percentsamples"), type="double", default=50, 
              help="Genes expressed with the given cutoff in atleast Nsamples should be considered 
                for this analysis. 
                Default is 50%"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output directory")
)

##Parse the command line arguments
opt_parser <- OptionParser(usage = "Rscript %prog [options] -f gene_expression.txt -d sample_datatable.txt 
       Rscript %prog -screen screen_name [options] -f gene_expression.txt -d sample_datatable.txt\n", 
  option_list=option_list, description = "  We use ridge regression to predict the target gene activity from the 
  binding TFs activity and then use a permutation test to assign significance to ridge regression coefficients. 
  The regularizer in the ridge regression will be chosen by 10-fold cross-validation. The target gene activity 
  will be randomly permuted followed by model fitting by ridge regression to obtain a null distribution of the 
  regression coefficients.
  Required inputs: 1) Gene expression file and 2) Datatable with sample information.
  Check help for detailed description of files")
opt <- parse_args(opt_parser)
gfile <- opt$gfile
dt.location <- opt$datatable
Nrandom <- opt$Nrandom
expressioncutoff <- opt$expressioncutoff
percentsamples <- opt$percentsamples
outdir <- opt$outdir

##Check for the correct input from commandline arguments
if (is.null(gfile)){
	print_help(opt_parser)
	stop("Gene expression file (RPKM/FPKM) must be supplied (input file)", call.=FALSE)
}

if (is.null(dt.location)){
	print_help(opt_parser)
	stop("Data table with sample information (input file)", call.=FALSE)
}

if (is.null(expressioncutoff)){
	expressioncutoff <- 1
}

if (is.null(percentsamples)){
	percentsamples <- 50
}

##Create the output directory or check for already present directory
if(!is.null(outdir)){
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }  
} else {
  outdir <- getwd()
}

##Check and Read the datatable
if(!file.exists(dt.location)){
  stop("Datatable file with sample information does not exist")
}
df.samples <- read.delim(dt.location, sep = "\t", header = T, stringsAsFactors = F)
sample.cols <- colnames(df.samples)
if(!("SAMPLE_NAME" %in% sample.cols)){
  stop("SAMPLE_NAME column missing in the datatable")
}
if(!("FILE" %in% sample.cols)){
  stop("SAMPLE_NAME column missing in the datatable")
}
expt.samples <- df.samples$SAMPLE_NAME
Nsamples <- length(expt.samples)
if(!(length(unique(expt.samples)) == Nsamples)){
  stop("Sample names in the datatable are not unique")
}

##Set the number of randomizations
##If no value for randomization is given the set it to 10000
if (is.null(Nrandom)){
  possible.perms <- factorial(Nsamples)
  if (possible.perms > 10000){
    possible.perms <- 10000
  }
}
##If the user has given a value of randomzations then makes sure that
##it is less than the number of possible randomizations
if (!is.null(Nrandom)){
  possible.perms <- factorial(Nsamples)
  if (possible.perms <- Nrandom){
    possible.perms <- Nrandom
  }
}

##Check and Read the gene expression file and filter all the samples present in the
##the datatable as well as the gene expression file
if(!file.exists(gfile)){
  stop("Genexpression file does not exist")
}
df.rpkm <- read.delim(gfile, sep = "\t", header = T, stringsAsFactors = F)
exprs.cols <- colnames(df.rpkm)
if(!("GENE_SYMBOL" %in% exprs.cols)){
  stop("GENE_SYMBOL column missing in the expression file")
}
exprs.samples <- setdiff(colnames(df.rpkm), "GENE_SYMBOL")
common.samples <- intersect(expt.samples, exprs.samples)
if(length(common.samples) < 2){
  stop("Needs atleast gene expression of two samples")
}

##Get all the expressed genes
exprs.genes <- toupper(getExpressedGenes(df.rpkm, samples = common.samples,
  present.prop = (percentsamples/100), rpkm.cutoff = expressioncutoff))

##Create a motif hit matrix
m.motif.hit <- getMotifHitMatrix(df.samples, samples = common.samples)
m.motif.hit <- m.motif.hit[intersect(exprs.genes, rownames(m.motif.hit)), ]

##All target genes
target.genes <- toupper(getALLTargetGenes(df.samples, samples = common.samples))

##Get the exprsd target genes and expressed tf
exprs.target.genes <- intersect(exprs.genes, target.genes)
exprs.tfs <- intersect(colnames(m.motif.hit), exprs.genes)

##Run ridge
ls.target.pvals <- sapply(exprs.target.genes, fitRidge, df.rpkm, 
  samples = common.samples, m.motif.hit = m.motif.hit, Nrand = possible.perms, simplify = F)
robjectfile <- file.path(outdir, "tf_target_ridge.Rds")
saveRDS(ls.target.pvals, file = robjectfile)

##Read the ridge results and get the significant TFs for each target gene
ls.target.pvals <- readRDS(robjectfile)
df.interactions <- getAllInteractions(ls.target.pvals)
interactionfile <- file.path(outdir, "tf_target_interactions.txt")
write.table(df.interactions, file = interactionfile, col.names = T, row.names = F, quote = F, sep = "\t")




















