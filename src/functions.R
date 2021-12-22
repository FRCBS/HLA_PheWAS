
## libs
sapply(c('tidyverse','data.table','ggpubr','viridis','rasterpdf', 'MASS', 'icd','mutoss', 'metap', 'cfdr.pleio',
         'sfsmisc',  'ghibli', 'RColorBrewer', 'SPAtest', 'ggstance', 'ggrepel', 'readxl', 'biomaRt' ,'qqplotr'),
       library, character.only=T)
# , 'Biostrings', 'circlize',
options(stringsAsFactors=F)


## --------------------------------------------------------
## functions
## --------------------------------------------------------

# beta difference Z test (two-way)
Z.test.2 <- function(m1, m2, se1, se2) { # means and sd's
  zstat <- (m1-m2) / sqrt( (se1^2) + (se2^2) )
  pnorm(abs(zstat), lower.tail=F)*2
}

# calculates protein divergence between two HLA alleles
HLAdivergence <- function(alleles, divergence.table, aligned.fasta, method='gms') {
  
  allele.match <- match(alleles, names(aligned.fasta))
  
  if (any(is.na(allele.match))) { NA } else {
    
    ind <- match(alleles, names(aligned.fasta))
    
    tmp <- data.frame(S1=str_split(toString(aligned.fasta[[ ind[1] ]]), '')[[1]],
                      S2=str_split(toString(aligned.fasta[[ ind[2] ]]), '')[[1]]) %>% unite(., 'AApair', S1, S2, sep='')
    tmp <- left_join(tmp, divergence.table, by='AApair') %>% na.omit
    
    if(method=='gms')      tmp <- tmp$gms_mean %>% sum
    if(method=='grantham') tmp <- tmp$grantham %>% sum
    if(method=='miyata')   tmp <- tmp$miyata %>% sum
    if(method=='sneath')   tmp <- tmp$sneath %>% sum
    
    tmp / length(aligned.fasta[[1]])
  }
}

# HLA notation formatting
formHLA <- function(x) {
  out <- unite(str_split_fixed(x, '\\.', 2) %>% data.frame, hla, X1, X2, sep='*')$hla
  out <- unite(str_split_fixed(out, '\\.', 2) %>% data.frame, hla, X1, X2, sep=':')$hla
  gsub(':$', '', out)
}


# from https://rdrr.io/github/alexploner/condFDR/src/R/condFDR.R
ccFDR = function(data, p1, p2, p_threshold = 1E-3, mc.cores = 1)
{
  ## Extract the data
  if ( !missing(data) ) {
    stopifnot( is.matrix(data) | is.data.frame(data)  )
    p1 = data[, p1]
    p2 = data[, p2]
  } else {
    stopifnot( length(p1) == length(p2))
    p1_name = deparse(substitute(p1))
    p2_name = deparse(substitute(p2))
  }
  
  ## Check: probabilities, no missing values
  doCheck = function(x) is.numeric(x) & !any(is.na(x)) & !any(x<0) & !any(x>1)
  stopifnot( doCheck(p1) & doCheck(p2) )
  
  ## Subset
  p_threshold = rep_len(p_threshold, length.out = 2)
  ndx = (p1 <= p_threshold[1]) | (p2 <= p_threshold[2])
  stopifnot( any(ndx) )
  p1  = p1[ndx]
  p2  = p2[ndx]
  
  ## Loop
  nn = length(p1)
  calc.denoms = function(i)
  {
    ## The edge point
    x = p1[i]
    y = p2[i]
    
    ## Vectors
    dd1 = p2 <= y
    dd2 = p1 <= x
    ee = dd1 & dd2
    
    ## Combine
    ee_n = length(which(ee))
    denom1 =  ee_n / length(which(dd2))
    denom2 =  ee_n / length(which(dd2))
    
    c(denom1, denom2)
  }
  ## This needs a bit care: unlist generates a vector with alternating denom1/denom2
  ## We can pour this vector into a matrix with two rows, without having to re-order,
  ## and can extract the rows in the next step
  denoms = matrix( unlist( parallel::mclapply(1:nn, calc.denoms, mc.cores = mc.cores) ), nrow = 2)
  
  ## Calibrate the p-values; note: use rows!
  cfdr1 = p1 / denoms[1, ]
  cfdr2 = p2 / denoms[2, ]
  ccfdr = pmax(cfdr1, cfdr2)
  
  ## Build output
  if ( !missing(data) ) {
    ret = cbind(data[ndx, ], cFDR1 = cfdr1, cFDR2 = cfdr2, ccFDR = ccfdr)
  } else {
    ret = data.frame(p1, p2, cFDR1 = cfdr1, cFDR2 = cfdr2, ccFDR = ccfdr)
    colnames(ret)[1:2] = c(p1_name, p2_name)
  }
  
  ret
}



