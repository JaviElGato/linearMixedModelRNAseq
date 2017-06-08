library(variancePartition)
library(doParallel)
# cl = makeCluster(8)
# registerDoParallel(cl)

foldChangeTest = function(Matrix, object, test, howMany, set1, set2, y, onWhich, setSubA, setSubB, colData, design, outplot, allowParallel=FALSE){
	# Matrix should be already normalized
	# object "edger" or "deseq2"
	# howMany groups are present in the dataset: 1 or 2
	# set1 and set2 corresponds to the number  INDIVIDUALS per group
	# test = FoldChanges Matrix; t-test or glm Regression ; variance Partition (Hoffman G. 2016 BMC Bioinformatics)
	# onWhich = On which group to apply the analysis: set1, set2, both
	# setSubA; setSubB = if "onWhich" is used, then need to specify the Subset number of samples
	# y = if asking to a GLM, the dependent variable should be specified
	# colData = variables to create the model - Used in Variance Partition but potentially to be used in linerMixedModels
	# design = Design model to be used in variancePartition or Mixed Models
	# outplot = name of Variance Explained Plot from variancePartition

	if(allowParallel){
		cl = makeCluster(8)
		registerDoParallel(cl)
	}

	if(object == "edger"){
		cat("object: edgeR \n")

		# Removing genes with CPM = 0 in at least one of the samples
		genesToKeep = rowSums(cpm(Matrix) > 1 ) >= 0.5 * ncol(Matrix)
		MatrixClean = DGEList(counts=Matrix[genesToKeep,])
		# TMM normalization
		MatrixClean = calcNormFactors(MatrixClean)
		# Voom normalization > Estimate precision weights for each gene and sample. This models uncertainty in expression measurments
		Matrix = voom(MatrixClean)
		Matrix = Matrix$E

	} else if (object == "deseq2") {
		cat("object: DESeq2 \n")

		# USE THIS PART FOR variancePartion
		# Estimate library size and correction scaling factors
		Matrix = estimateSizeFactors(Matrix)
		# Removing genes with FPKM = 0 in at least one of the samples
		genesToKeep = rowSums(fpm(Matrix) > 1 ) >= 1 * ncol(Matrix)
		# Compute Log2 Fragments per million
		Matrix = log2(fpm(Matrix)[genesToKeep,] + 1)

	}

	numberGenes = dim(Matrix)[1]
	numberSamples = dim(Matrix)[2]
	cat(numberGenes, "Genes in", numberSamples, "samples", "\n")

############################################################
	
	returnFoldChangesPerIndv = function(matrixIn, numberGenes, numberSamples, set1, set2, howMany){
		

		if (howMany > 1) { 
			
			resultsMatrix = matrix(seq(1,numberGenes),  nrow = numberGenes, ncol = set1+set2)
			rownames(resultsMatrix) = rownames(matrixIn)
			for (col in 1:set1){
			# print(col)
				for(gene in 1:numberGenes){
				# cat(col, gene, "\n")
					valuePre = matrixIn[gene, col]
					valuePost = matrixIn[gene, (col+set1)]
					# cat(valuePre, valuePost, "\n")
					# foldChange = abs((valuePost/valuePre)-1)
					foldChange = log2(valuePost/valuePre)
					resultsMatrix[gene,col] = foldChange
				}

			}

			for (col in ((set1*2)+1):(numberSamples-set2)) {
				for(gene in 1:numberGenes){
					# cat(col, gene, "\n")
					valuePre = matrixIn[gene, col]
					valuePost = matrixIn[gene, (col+set2)]
					# cat(valuePre, valuePost, "\n")
					# foldChange = abs((valuePost/valuePre)-1)
					foldChange = log2(valuePost/valuePre)
					resultsMatrix[gene,col-set1] = foldChange
				}
			
			}

		} else {

			resultsMatrix = matrix(seq(1,numberGenes),  nrow = numberGenes, ncol = set1)
			rownames(resultsMatrix) = rownames(matrixIn)

			for (col in 1:set1){
			# print(col)
				for(gene in 1:numberGenes){
				# cat(col, gene, "\n")
					valuePre = matrixIn[gene, col]
					valuePost = matrixIn[gene, (col+set1)]
					# cat(valuePre, valuePost, "\n")
					# foldChange = abs((valuePost/valuePre)-1)
					foldChange = log2(valuePost/valuePre)
					resultsMatrix[gene,col] = foldChange
				}

			}

		} 

		resultsMatrix

	}

############################################################

############################################################

	returnTtest = function(foldChangesMatrix, set1){

		numberGenes = dim(foldChangesMatrix)[1]
		numberIndiv = dim(foldChangesMatrix)[2]
		t_testResults = data.frame(mean1=rep(1,numberGenes), mean2=rep(1,numberGenes), CI95_L=rep(1,numberGenes), CI95_U=rep(1,numberGenes), t=rep(1,numberGenes), df=rep(1,numberGenes), pvalue=rep(1,numberGenes), pvalue.adj=rep(1,numberGenes))
		rownames(t_testResults) = rownames(foldChangesMatrix)

		for(gene in 1:numberGenes){

			toprint = gene %% 1000
			if(toprint == 0){
				cat(paste(gene, ",", sep=""))
			}

			geneTtestResult = t.test(foldChangesMatrix[gene,1:set1], foldChangesMatrix[gene,(set1+1):numberIndiv] )
			mean1 = geneTtestResult$estimate[[1]]
			mean2 = geneTtestResult$estimate[[2]]
			ci95L = geneTtestResult$conf.int[[1]]
			ci95U = geneTtestResult$conf.int[[2]]
			t = geneTtestResult$statistic[[1]]
			df = geneTtestResult$parameter[[1]]
			pval = geneTtestResult$p.value

			t_testResults[gene,] = c(mean1=mean1, mean2=mean2, CI95_L=ci95L, CI95_U=ci95U, t=t, df=df, pvalue=pval, pvalue.adj=0 )

		}

		cat("Done! \n")

		pvaladj = p.adjust(t_testResults$pvalue, , method="BH")
		t_testResults[,8] =  pvaladj
		
		t_testResults



	}

############################################################


############################################################
	
	returnGlm = function(foldChangesMatrix, y, set1, set2){

		numberGenes = dim(foldChangesMatrix)[1]
		numberIndiv = dim(foldChangesMatrix)[2]
		glm_testResults = data.frame(beta=rep(1,numberGenes), CI95_L=rep(1,numberGenes), CI95_U=rep(1,numberGenes), pvalue=rep(1,numberGenes), pvalue.adj=rep(1,numberGenes))
		rownames(glm_testResults) = rownames(foldChangesMatrix)

		for(gene in 1:numberGenes){

			toprint = gene %% 1000
			if(toprint == 0){
				cat(paste(gene, ",", sep=""))
			}

			geneGlm = glm(foldChangesMatrix[gene,] ~ y )
			geneGlmConfint = suppressMessages(confint(geneGlm))
			geneGlmSumm = summary(glm(y ~ foldChangesMatrix[gene,]))
			beta = geneGlmSumm$coefficients[[2]]
			ci95L = geneGlmConfint[[2]]
			ci95U = geneGlmConfint[[4]]
			pvalue = geneGlmSumm$coefficients[[8]]

			glm_testResults[gene,] = c(beta=beta, CI95L=ci95L, CI95U=ci95U, pvalue=pvalue, pvalue.adj=0 )

		}

		cat("Done! \n")

		pvaladj = p.adjust(glm_testResults$pvalue, method="BH")
		glm_testResults[,5] =  pvaladj
		
		glm_testResults
	}

############################################################

############################################################

	returnVariancePartition = function(matrixIn, colData, design, outplot){

		# Fit the Model
		varPart = fitExtractVarPartModel(exprObj=matrixIn, data=colData, formula=design)

		# Sort cols
		varPartSorted = sortCols(varPart)

		# Plot Variance Explained
		varianceExplainedPlot = plotVarPart(varPartSorted)
		ggsave(outplot, varianceExplainedPlot)

		varPartSorted

	}
############################################################

	if (test == "FoldChanges"){

		returnFoldChangesPerIndv(matrixIn=Matrix, numberGenes=numberGenes, numberSamples=numberSamples, set1=set1, set2=set2, howMany=howMany)

	} else if (test == "t-test"){
		cat("Calculating t-test... \n")
		if(howMany > 1){
			foldChangesMatrix = returnFoldChangesPerIndv(matrixIn=Matrix, numberGenes=numberGenes, numberSamples=numberSamples, set1=set1, set2=set2)
			returnTtest(foldChangesMatrix=foldChangesMatrix , set1=set1)
		} else {
			returnTtest(foldChangesMatrix=Matrix , set1=set1)
		}
		

	} else if (test == "glm"){
		cat("Calculating GLM... \n")
		if(howMany > 1){
			foldChangesMatrix = returnFoldChangesPerIndv(matrixIn=Matrix, numberGenes=numberGenes, numberSamples=numberSamples, set1=set1, set2=set2)
			returnGlm(foldChangesMatrix, y, set1, set2)
		} else {
			returnGlm(Matrix, y, set1, set2)

		}


	} else if (test == "variancePartition"){ ####### CHECK AND FIX
		cat(paste("variancePartition analysis with design:", design, sep=" "), "\n")

		if ( onWhich == "both" ){
			
			cat("variancePartition in the full matrix \n")
			returnVariancePartition(matrixIn=Matrix, colData=colData, design=design, outplot=outplot)	
			
		} else if (onWhich == "set1") {

			cat("variancePartition only on Set1 \n")
			MatrixSubset = Matrix[,1:(set1*2)]
			returnVariancePartition(matrixIn=MatrixSubset, colData=colData, design=design, outplot=outplot)	
			# cat("Under development \n")

		} else if ((onWhich == "set2")) {

			cat("variancePartition only on Set2 \n")
			numberIndiv = dim(Matrix)[2]
			MatrixSubset = Matrix[,((set1*2)+1):numberIndiv]
			print(head(MatrixSubset))
			# MatrixSubset = Matrix[,12:32]
			cat(numberIndiv, dim(MatrixSubset), "\n")
			returnVariancePartition(matrixIn=MatrixSubset, colData=colData, design=design, outplot=outplot)	
			# cat("Under development \n")

		}


	}



}
