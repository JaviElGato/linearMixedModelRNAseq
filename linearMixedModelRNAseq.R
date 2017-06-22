library(variancePartition)
library(doParallel)
library(lme4)
# cl = makeCluster(8)
# registerDoParallel(cl)

foldChangeTest = function(Matrix, object, test, howMany, set1, set2, y, onWhich, setSubA, setSubB, colData, design, outplot, allowParallel=FALSE, designFull, designNull){
	# Matrix should be already normalized
	# object "edger", "deseq2", "normCPM"
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
		cl = makeCluster(4)
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

	} else if ( object == "normCPM"){
		cat("object: Matrix of normalized CPM\n")

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
	
	returnLm = function(foldChangesMatrix, y, set1, set2, designNull, designFull, colData){

		numberGenes = dim(foldChangesMatrix)[1]
		numberIndiv = dim(foldChangesMatrix)[2]
		# glm_testResults = data.frame(beta=rep(1,numberGenes), CI95_L=rep(1,numberGenes), CI95_U=rep(1,numberGenes), varExplained=rep(1,numberGenes), pvalue=rep(1,numberGenes), pvalue.adj=rep(1,numberGenes))
		glm_testResults = data.frame(beta=rep(1,numberGenes), CI95_L=rep(1,numberGenes), CI95_U=rep(1,numberGenes), varExplainedNull=rep(1,numberGenes), varExplainedFull=rep(1,numberGenes), nullAIC=rep(1,numberGenes), fullAIC=rep(1,numberGenes), diffAIC=rep(1,numberGenes), pvalue=rep(1,numberGenes), pvalue.adj=rep(1,numberGenes))
		rownames(glm_testResults) = rownames(foldChangesMatrix)

		lmp = function (modelobject) {
			   if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
			   f <- summary(modelobject)$fstatistic
			   p <- pf(f[1],f[2],f[3],lower.tail=F)
			   attributes(p) <- NULL
			   return(p)
			}

		for(gene in 1:numberGenes){

			toprint = gene %% 1000
			if(toprint == 0){
				cat(paste(gene, ",", sep=""))
			}

			# designFormula = as.formula(paste("foldChangesMatrix[gene,] ~", design))
			designFormulaNull = as.formula(paste("foldChangesMatrix[gene,] ~", designNull))
			designFormulaFull = as.formula(paste("foldChangesMatrix[gene,] ~", designFull))
			# print(designFormula)
			geneGlmNull = lm(designFormulaNull, data=colData)
			geneGlmFull = lm(designFormulaFull, data=colData)
			# geneGlm = glm(designFormula, data=colData, family="gaussian" )
			# geneGlm = glm(foldChangesMatrix[gene,] ~ design, data=colData, family="gaussian" )
			geneGlmConfint = suppressMessages(confint(geneGlmFull))

			geneGlmNullSumm = summary(geneGlmNull)
			geneGlmFullSumm = summary(geneGlmFull)
			# geneGlmSumm = summary(glm(designFormula, data=colData, family="gaussian" ))
			# geneGlmSumm = summary(glm(foldChangesMatrix[gene,] ~ design, data=colData, family="gaussian" ))
			coefficients = dim(geneGlmFullSumm$coefficients)
			# coefficients = length(geneGlmSumm$coefficients[,1])
			# beta = geneGlmSumm$coefficients[[2]]
			beta = geneGlmFullSumm$coefficients[,1][[coefficients[1]]]
			ci95L = geneGlmConfint[[2]]
			ci95U = geneGlmConfint[[4]]
			varianceExplainedNull = geneGlmNullSumm$r.squared # $adj.r.squared
			varianceExplainedFull = geneGlmFullSumm$r.squared

			null_AIC = AIC(geneGlmNull)
			full_AIC = AIC(geneGlmFull)
			diff_AIC = abs(null_AIC - full_AIC)

			pvalue = lmp(geneGlmFull)

			# cat(beta, ci95L, ci95U, null_AIC, full_AIC, diff_AIC, pvalue, "\n\n")
			# glm_testResults[gene,] = c(beta=beta, CI95L=ci95L, CI95U=ci95U, varExplained=varianceExplained, pvalue=pvalue, pvalue.adj=0 )
			glm_testResults[gene,] = c(beta=beta, CI95L=ci95L, CI95U=ci95U, varExplainedNull=varianceExplainedNull, varExplainedFull=varianceExplainedFull, nullAIC=null_AIC, fullAIC=full_AIC, diffAIC=diff_AIC, pvalue=pvalue, pvalue.adj=0 )

		}

		cat("Done! \n")

		pvaladj = p.adjust(glm_testResults$pvalue, method="BH")
		glm_testResults[,10] =  pvaladj
		
		glm_testResults
	}

############################################################

	returnLinearMixedModel = function(foldChangesMatrix, designNull, designFull, colData, ciMethod){

		numberGenes = dim(foldChangesMatrix)[1]
		numberIndiv = dim(foldChangesMatrix)[2]
		lmm_testResults = data.frame(beta=rep(1,numberGenes), CI95_L=rep(1,numberGenes), CI95_U=rep(1,numberGenes), nullAIC=rep(1,numberGenes), fullAIC=rep(1,numberGenes), diffAIC=rep(1,numberGenes), pvalue=rep(1,numberGenes), pvalue.adj=rep(1,numberGenes))
		rownames(lmm_testResults) = rownames(foldChangesMatrix)


		for(gene in 1:numberGenes){
			# cat(paste(gene, ",", sep=""))
			toprint = gene %% 100
			if(toprint == 0){
				cat(paste(gene, ",", sep=""))
			}

			designFormulaNull = as.formula(paste("foldChangesMatrix[gene,] ~", designNull))
			designFormulaFull = as.formula(paste("foldChangesMatrix[gene,] ~", designFull))
			# print(designFormula)
			geneLmmNull = lmer(designFormulaNull, data=colData)
			geneLmmFull = lmer(designFormulaFull, data=colData)
			# geneGlm = glm(designFormula, data=colData, family="gaussian" )
			# geneGlm = glm(foldChangesMatrix[gene,] ~ design, data=colData, family="gaussian" )
			geneLmmConfint = suppressMessages(confint(geneLmmFull, method=ciMethod))
			cis = length(geneLmmConfint[,1])

			geneLmmSumm = summary(geneLmmFull)
			# geneGlmSumm = summary(glm(designFormula, data=colData, family="gaussian" ))
			# geneGlmSumm = summary(glm(foldChangesMatrix[gene,] ~ design, data=colData, family="gaussian" ))
			coefficients = length(geneLmmSumm$coefficients[,1])
			# beta = geneGlmSumm$coefficients[[2]]
			beta = geneLmmSumm$coefficients[,1][[coefficients]]

			ci95L = geneLmmConfint[,1][[cis]]
			ci95U = geneLmmConfint[,2][[cis]]
			
			# pvalue = geneGlmSumm$coefficients[,4][[coefficients]]
			geneLmmAnova = suppressMessages(anova(geneLmmNull, geneLmmFull))
			null_AIC = geneLmmAnova$AIC[1]
			full_AIC = geneLmmAnova$AIC[2]
			diff_AIC = abs(null_AIC - full_AIC)
			pvalue = geneLmmAnova[[8]][2]
			
			
			lmm_testResults[gene,] = c(beta=beta, CI95L=ci95L, CI95U=ci95U, nullAIC=null_AIC, fullAIC=full_AIC, diffAIC=diff_AIC, pvalue=pvalue, pvalue.adj=0 )
		

		}

		cat("Done! \n")

		pvaladj = p.adjust(lmm_testResults$pvalue, method="BH")
		lmm_testResults[,8] =  pvaladj
		
		lmm_testResults
	}


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
			foldChangesMatrix = returnFoldChangesPerIndv(matrixIn=Matrix, numberGenes=numberGenes, numberSamples=numberSamples, set1=set1, set2=set2, howMany=howMany)
			returnTtest(foldChangesMatrix=foldChangesMatrix , set1=set1)
		} else {
			returnTtest(foldChangesMatrix=Matrix , set1=set1)
		}
		

	} else if (test == "lm"){
		cat("Calculating LM... \n")
		if(howMany > 1){
			foldChangesMatrix = returnFoldChangesPerIndv(matrixIn=Matrix, numberGenes=numberGenes, numberSamples=numberSamples, set1=set1, set2=set2, howMany=howMany)
			returnLm(foldChangesMatrix, y, set1, set2, design, colData=colData)
		} else {
			returnLm(Matrix, y, set1, set2,  designNull, designFull, colData=colData)

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


	} else if (test == "lmm"){
		cat("Calculating LMM... \n")

		if(howMany > 1){
			foldChangesMatrix = returnFoldChangesPerIndv(matrixIn=Matrix, numberGenes=numberGenes, numberSamples=numberSamples, set1=set1, set2=set2, howMany=howMany)
			returnLinearMixedModel(foldChangesMatrix, designNull, designFull, colData, ciMethod="Wald")
		} else {
			returnLinearMixedModel(Matrix, designNull, designFull, colData, ciMethod="Wald")
		}

		
	}



}
