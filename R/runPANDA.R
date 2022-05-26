runPANDA <- function(motif,expr=NULL,ppi=NULL,alpha=0.1,hamming=0.001,
                  iter=NA,output=c('regulatory','coexpression','cooperative'),
                  zScale=TRUE,progress=TRUE,randomize=c("None", "within.gene", "by.gene"), assoc.method="pearson",
                  scale.by.present=FALSE,edgelist=FALSE,remove.missing.ppi=FALSE,
                  remove.missing.motif=FALSE,remove.missing.genes=FALSE,mode="intersection"){

  randomize <- match.arg(randomize)
  if(progress)
    print('Initializing and validating')
#
#   if(class(expr)=="ExpressionSet")
#     expr <- assayData(expr)[["exprs"]]

  if (is.null(expr)){
    # Use only the motif data here for the gene list
    num.conditions <- 0
    if (randomize!="None"){
      warning("Randomization ignored because gene expression is not used.")
      randomize <- "None"
    }
  } else {
    if(mode=='legacy'){
      if(remove.missing.genes){
        # remove genes from expression data that are not in the motif data
        n <- nrow(expr)
        expr <- expr[which(rownames(expr)%in%motif[,2]),]
        message(sprintf("%s genes removed that were not present in motif", n-nrow(expr)))
      }
      if(remove.missing.motif){
        # remove genes from motif data that are not in the expression data
        n <- nrow(motif)
        motif <- motif[which(motif[,2]%in%rownames(expr)),]
        message(sprintf("%s motif edges removed that targeted genes missing in expression data", n-nrow(motif)))
      }
      # Use the motif data AND the expr data (if provided) for the gene list
      # Keep everything sorted alphabetically
      expr <- expr[order(rownames(expr)),]
    } else if(mode=='union'){
      gene.names=unique(union(rownames(expr),unique(motif[,2])))
      tf.names  =unique(union(unique(ppi[,1]),unique(motif[,1])))
      num.TFs    <- length(tf.names)
      num.genes  <- length(gene.names)
      # gene expression matrix
      expr1=as.data.frame(matrix(0,num.genes,ncol(expr)))
      rownames(expr1)=gene.names
      expr1[which(gene.names%in%rownames(expr)),]=expr[]
      expr=expr1
      #PPI matrix
      tfCoopNetwork <- matrix(0,num.TFs,num.TFs)
      colnames(tfCoopNetwork)=tf.names
      rownames(tfCoopNetwork)=tf.names
      Idx1 <- match(ppi[,1], tf.names);
      Idx2 <- match(ppi[,2], tf.names);
      Idx <- (Idx2-1)*num.TFs+Idx1;
      tfCoopNetwork[Idx] <- ppi[,3];
      Idx <- (Idx1-1)*num.TFs+Idx2;
      tfCoopNetwork[Idx] <- ppi[,3];
      #Motif matrix
      regulatoryNetwork=matrix(0,num.TFs,num.genes)
      colnames(regulatoryNetwork) = gene.names
      rownames(regulatoryNetwork) = tf.names
      Idx1=match(motif[,1], tf.names);
      Idx2=match(motif[,2], gene.names);
      Idx=(Idx2-1)*num.TFs+Idx1;
      regulatoryNetwork[Idx]=motif[,3]
    } else if(mode=='intersection'){

      gene.names = sort(intersect(rownames(expr), motif[,2]))
      tf.names  = sort(intersect(unique(c(ppi[,1], ppi[,2])), motif[,1]))
      num.TFs   = length(tf.names)
      num.genes  = length(gene.names)

      # Gene expression matrix
      expr = expr[gene.names,]

      # PPI matrix
      ppi <- ppi[ppi[,1] %in% tf.names,]
      ppi <- ppi[ppi[,2] %in% tf.names,]
      tfCoopNetwork <- Matrix::sparseMatrix(i = as.numeric(factor(ppi[,1], tf.names)),
                                            j = as.numeric(factor(ppi[,2], tf.names)),
                                            x = ppi[,3],
                                            dims = c(num.TFs, num.TFs),
                                            dimnames = list(tf.names, tf.names))
      tfCoopNetwork <- as.matrix(tfCoopNetwork)


      #Motif matrix
      motif <- motif[motif[,1] %in% tf.names,]
      motif <- motif[motif[,2] %in% gene.names,]
      regulatoryNetwork <- Matrix::sparseMatrix(i = as.numeric(factor(motif[,1], tf.names)),
                                               j = as.numeric(factor(motif[,2], gene.names)),
                                               x = motif[,3],
                                               dims = c(num.TFs, num.genes),
                                               dimnames = list(tf.names, gene.names))
      regulatoryNetwork <- as.matrix(regulatoryNetwork)
    }

    num.conditions <- ncol(expr)
    if (randomize=='within.gene'){
      expr <- t(apply(expr, 1, sample))
      if(progress)
        print("Randomizing by reordering each gene's expression")
    } else if (randomize=='by.gene'){
      rownames(expr) <- sample(rownames(expr))
      expr           <- expr[order(rownames(expr)),]
      if(progress)
        print("Randomizing by reordering each gene labels")
    }
  }

  if (mode=='legacy'){
    # Create vectors for TF names and Gene names from motif dataset
    tf.names   <- sort(unique(motif[,1]))
    gene.names <- sort(unique(rownames(expr)))
    num.TFs    <- length(tf.names)
    num.genes  <- length(gene.names)
  }

  # Bad data checking
  if (num.genes==0){
    stop("Error validating data.  No matched genes.\n  Please ensure that gene names in expression data match gene names in motif data")
  }

  if(num.conditions==0) {
    warning('No expression data given.  PANDA will run based on an identity co-regulation matrix')
    geneCoreg <- diag(num.genes)
  } else if(num.conditions<3) {
    warning('Not enough expression conditions detected to calculate correlation. Co-regulation network will be initialized to an identity matrix.')
    geneCoreg <- diag(num.genes)
  } else {

    if(scale.by.present){
      num.positive=(expr>0)%*%t((expr>0))
      if(assoc.method == 'pcr'){
        geneCoreg <- as.matrix(pcNet(expr)) * (num.positive/num.conditions)
      } else {
        geneCoreg <- cor(t(expr), method=assoc.method, use="pairwise.complete.obs")
      }

    } else {
      if(assoc.method == 'pcr'){
        geneCoreg <- as.matrix(pcNet(expr))
      } else {
        geneCoreg <- cor(t(expr), method=assoc.method, use="pairwise.complete.obs")
      }
    }
    if(progress)
      print('Verified sufficient samples')
  }
  if (any(is.na(geneCoreg))){ #check for NA and replace them by zero
    diag(geneCoreg)=1
    geneCoreg[is.na(geneCoreg)]=0
  }

  if (any(duplicated(motif))) {
    warning("Duplicate edges have been found in the motif data. Weights will be summed.")
    motif <- aggregate(motif[,3], by=list(motif[,1], motif[,2]), FUN=sum)
  }

  # Prior Regulatory Network
  if(mode=='legacy'){
    Idx1=match(motif[,1], tf.names);
    Idx2=match(motif[,2], gene.names);
    Idx=(Idx2-1)*num.TFs+Idx1;
    regulatoryNetwork=matrix(data=0, num.TFs, num.genes);
    regulatoryNetwork[Idx]=motif[,3]
    colnames(regulatoryNetwork) <- gene.names
    rownames(regulatoryNetwork) <- tf.names
    # PPI data
    # If no ppi data is given, we use the identity matrix
    tfCoopNetwork <- diag(num.TFs)
    # Else we convert our two-column data.frame to a matrix
    if (!is.null(ppi)){
      if(any(duplicated(ppi))){
        warning("Duplicate edges have been found in the PPI data. Weights will be summed.")
        ppi <- aggregate(ppi[,3], by=list(ppi[,1], ppi[,2]), FUN=sum)
      }
      if(remove.missing.ppi){
        # remove edges in the PPI data that target TFs not in the motif
        n <- nrow(ppi)
        ppi <- ppi[which(ppi[,1]%in%tf.names & ppi[,2]%in%tf.names),]
        message(sprintf("%s PPI edges removed that were not present in motif", n-nrow(ppi)))
      }
      Idx1 <- match(ppi[,1], tf.names);
      Idx2 <- match(ppi[,2], tf.names);
      Idx <- (Idx2-1)*num.TFs+Idx1;
      tfCoopNetwork[Idx] <- ppi[,3];
      Idx <- (Idx1-1)*num.TFs+Idx2;
      tfCoopNetwork[Idx] <- ppi[,3];
    }
    colnames(tfCoopNetwork) <- tf.names
    rownames(tfCoopNetwork) <- tf.names
  }

  ## Run PANDA ##
  tic=proc.time()[3]

  if(progress)
    print('Normalizing networks...')
  regulatoryNetwork = normalizeNetwork(regulatoryNetwork)
  tfCoopNetwork     = normalizeNetwork(tfCoopNetwork)
  geneCoreg         = normalizeNetwork(geneCoreg)

  if(progress)
    print('Learning Network...')

  minusAlpha = 1-alpha
  step=0
  hamming_cur=1
  if(progress)
    print("Using tanimoto similarity")
  while(hamming_cur>hamming){
    if ((!is.na(iter))&&step>=iter){
      print(paste("Reached maximum iterations, iter =",iter),sep="")
      break
    }
    Responsibility=tanimoto(tfCoopNetwork, regulatoryNetwork)
    Availability=tanimoto(regulatoryNetwork, geneCoreg)
    RA = 0.5*(Responsibility+Availability)

    hamming_cur=sum(abs(regulatoryNetwork-RA))/(num.TFs*num.genes)
    regulatoryNetwork=minusAlpha*regulatoryNetwork + alpha*RA

    ppi=tanimoto(regulatoryNetwork, t(regulatoryNetwork))
    ppi=update.diagonal(ppi, num.TFs, alpha, step)
    tfCoopNetwork=minusAlpha*tfCoopNetwork + alpha*ppi

    CoReg2=tanimoto(t(regulatoryNetwork), regulatoryNetwork)
    CoReg2=update.diagonal(CoReg2, num.genes, alpha, step)
    geneCoreg=minusAlpha*geneCoreg + alpha*CoReg2

    if(progress)
      message("Iteration", step,": hamming distance =", round(hamming_cur,5))
    step=step+1
  }

  toc=proc.time()[3] - tic
  if(progress)
    message("Successfully ran PANDA on ", num.genes, " Genes and ", num.TFs, " TFs.\nTime elapsed:", round(toc,2), "seconds.")
  prepResult(zScale, output, regulatoryNetwork, geneCoreg, tfCoopNetwork, edgelist, motif)
}
