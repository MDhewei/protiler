## The functions to filter signals of inefficient sgRNAs and Do segmentation on tiling screen signals
## @Date: 1/3/2019
## @authors: Han Xu, Wei He
## @email: whe3@mdanderson.org


library(breakfast) 
library(stringr)

args = commandArgs(trailingOnly=TRUE)  ## Use aruguments from user input with command line


# This function is to supress outliers of the signal

# Input: (1) x: the one-dimension vector of tiling screen scores
#        (2) win.size: the size of sliding window to detect outliers 
#        (3) th: the threshold used to determine outliers, defaut is 2.0

# Output: This function will return new list of signals with outlier ajusted to normal.

supress.outlier <- function(x,win.size=11, th=2.0)
{
  y <- runmed(x,win.size)  ## Compute running medians of odd span
  s <- mad(diff(x))/sqrt(2) ## Compute the median absolute deviation of lagged differences
  
  ## Adjust score to normal if the value two times divation larger or smaller than running median   
  z <- ifelse(x>y+th*s,y+th*s,ifelse(x<y-th*s,y-th*s,x)) 
    
  return(z)  ## Return adjusted signals
}


# This function is to check if an sgRNA is effecient by comparison with neighboring sgRNAs

# Input: (1) x: a numeric vector recording score of current sgRNA and its neribors
#        (2) quant: the probablity to produce sample quantiles

# Output: This function returns 0(inefficient) or 1(efficient)

is.efficient <- function(x, quant = 2/3)
{
  n <- (length(x)-1)/2
  #if ((x[n+1]<=median(x[1:n])) | (x[n+1]<=median(x[(n+2):length(x)])))
    
  ## quantile funciton produces sample quantiles corresponding to the given probabilities
  if ((x[n+1]<=quantile(x[1:n],quant)) | (x[n+1]<=quantile(x[(n+2):length(x)], quant)))
  {
    return(1)  ## The sgRNA is efficient if score smaller than 2/3 quantile value of neighbors in either side
  }
  else
  {
    return(0) ## The sgRNA is efficient if score larger than 2/3 quantile value of neighbors in both sides
  }
}


# This function segmentation of protein based on sgRNA values using TGUH method

# Input: (1) gene: the symbol of target gene. eg: CREBBP
#        (2) value: a vector recording functional scores for all the tiling sgRNAs
#        (3) AA.loc: numeric vector recording residue postions which sgRNAs cut
#        (4) half.size: the number of neighboring sgRNAs selected to filter inefficient guides
#        (5) th.outlier: the threshold used to determine outliers, defaut is 2.0
#        (6) th.tguh: the threshold used in TGUH method to detect changing points in input signals.

# Output: This function returns the following objects:
#         (1) v.filtered: the one-demension signals after filtering the outliers
#         (2) is.efficient: vector recording the efficiency of all the sgRNAs
#         (3) sgRNA.est: the segmented scores for all the sgRNAs.
#         (4) AA.est: the segmented scores along the amino acid positions
#         (5) segments: A table contains all the information for the segments called by TGUH.

protile.segmentation <- function(gene,value, AA.loc, half.size = 5, th.outlier = 2, th.tguh = 1.5)
{
  x <- value[order(AA.loc, decreasing = F)]  ##sort the input signals based on amino acids order(small to large)
  AA <- AA.loc[order(AA.loc,decreasing = F)] ##sort residue positions based on amino acids order(small to large)
  
  y <- c(rep(median(x),half.size),x,rep(median(x),half.size))
  
  ## Judge efficiency for all the sgRNAs
  flags <- sapply((half.size+1):(half.size+length(x)), function(i) is.efficient(y[(i-half.size):(i+half.size)]))
  
  ## Record all the efficient sgRNA index.
  index <- c() 
  
  for (i in 1:length(x))
  {
    if (flags[i])
    {
      index <- c(index,i)
    }
  }
      
  ## Supress outliers for those efficient sgRNAs
  w <- supress.outlier(x[which(flags==1)], th = th.outlier)
      
  ## Do Segmentation on filtered data using TGUH method.    
  seg <- segment.mean(w,th.const=th.tguh) 
  
  cpt <- c()
      
  ## Annotate changing point in the data
  if (length(seg$cpt) >0)
  {
    for (i in 1:length(seg$cpt))
    {
      if ((seg$cpt[i]==1)||(seg$cpt[i]==length(w)))
      {
        next
      }
      cpt <- c(cpt, floor((index[seg$cpt[i]]+index[seg$cpt[i]+1])/2))
    }
  }
  
  v.filtered <- rep(NA, length(x))
  
  for (i in 1:length(w))
  {
    v.filtered[index[i]] <- w[i]
  }
  
  sgRNA.est <- rep(mean(w),length(x))  ##initialize sgRNA score as means of input score
  AA.filtered <- AA[which(flags==1)]
  
  segments <- c()
  
  if (length(cpt) > 0)
  {
    m <-mean(v.filtered[1:cpt[1]],na.rm=T)
    sgRNA.est[1:cpt[1]] <- m   # assign the mean score of segments to different sgRNAs 
    
    segments <- rbind(segments, data.frame(start=1,
                                         end=cpt[1],
                                         AA.start=1,
                                         AA.end=floor((AA.filtered[seg$cpt[1]]+AA.filtered[seg$cpt[1]+1])/2),
                                         n=seg$cpt[1],
                                         m=m))

    
    if (length(cpt) > 1)  ## If we got multiple changing point
    {
      for (i in 1:(length(cpt)-1))
      {
        m <- mean(v.filtered[(cpt[i]+1):cpt[i+1]], na.rm =T)
        sgRNA.est[(cpt[i]+1):cpt[i+1]] <- m
        
        segments <- rbind(segments, data.frame(start=cpt[i]+1,
                                             end=cpt[i+1],
                                             AA.start=ceiling((AA.filtered[seg$cpt[i]]+AA.filtered[seg$cpt[i]+1])/2),
                                             AA.end=floor((AA.filtered[seg$cpt[i+1]]+AA.filtered[seg$cpt[i+1]+1])/2),
                                             n=seg$cpt[i+1]-seg$cpt[i],
                                             m=m))
      }
    }
    
    ## Get mean score of each segments
    m <- mean(v.filtered[(cpt[length(cpt)]+1):length(x)], na.rm=T)
    sgRNA.est[(cpt[length(cpt)]+1):length(x)] <- m
    
      
    ## Generate the table containing all the information for all the segments.
    segments <- rbind(segments, data.frame(start=cpt[length(cpt)]+1,end=length(x),
                                   AA.start=ceiling((AA.filtered[seg$cpt[length(cpt)]]+AA.filtered[seg$cpt[length(cpt)]+1])/2),
                                   AA.end=AA[length(AA)],n=length(w)-seg$cpt[length(cpt)],m=m))
   }
  else  ## If we got no changing points
  {
    m <- mean(v.filtered, na.rm=T)
    sgRNA.est[1:length(x)] <- m
    segments <- rbind(segments, data.frame(start=1,
                                           end=length(x),
                                           AA.start=1,
                                           AA.end=AA[length(AA)],
                                           n=length(w),
                                           m=m))
  }
  
  segments <- segments[order(segments$m),]
  
  ## Judge whether a segment is HS region
  is.HS.site <- rep(F, nrow(segments))
  is.HS.site[1] <- T
  
  if (nrow(segments)>2)
  for (i in 2:(nrow(segments)-1))
  {
    if (segments$m[i] > mean(v.filtered, na.rm=T))
    {
      break
    }
    
    m.lower <- weighted.mean(segments$m[1:(i-1)],segments$n[1:(i-1)])
    m.higher <- weighted.mean(segments$m[(i+1):nrow(segments)],segments$n[(i+1):nrow(segments)])
    
    if (segments$m[i]-m.lower < m.higher-segments$m[i])
    {
      is.HS.site[i] <- T
    }
    else
    {
        break
      }
  }
  
  segments <- cbind(segments, is.HS.site)
  segments <- segments[order(segments$AA.start), ]
  
  AA.est <- rep("NA",AA[length(AA)])

  for (i in 1:nrow(segments))
  {
    AA.est[segments$AA.start[i]:segments$AA.end[i]] <- segments$m[i]
  }
  
  return(list(v.filtered=v.filtered, is.efficient = flags, sgRNA.est=sgRNA.est, AA.est=AA.est, segments=segments))
}


# This function is to call the hyper sensitive regions in a protein encoded by certain gene

# Input: (1) inputfile: The table file recording the CRISPR tiling screen data
#        (2) gene: the symbol of target gene, eg:CREBBP
#        (3) col: the column number(s) of crispr scores 
#        (4) size: the number of neighboring sgRNAs selected to filter inefficient guides.
#        (5) th1: the threshold for surpress the outliers
#        (6) th2: the threshold for TGUH segmentation.

# Output: This function generate three table files for next step analysis and visualization.
#         (1) 'XX_Score.csv': Record the functional scores at each residue after filter
#         (2) 'XX_Est.csv': Record all the mean score for all the segments along the protein
#         (3) 'XX_Segments.csv': the segments call by TGUH together with other information for each segment

CallHSRegion <- function(inputfile,gene,col,size,th1,th2)
  {
   if (str_detect(inputfile,'.txt')==TRUE)
      {
      df.input <- read.table(inputfile,header=TRUE,sep="\t")
       }
   if (str_detect(inputfile,'.csv')==TRUE)
      {
      df.input <- read.csv(inputfile,header=TRUE,sep=",")
       }
   
   df.child <- df.input[df.input$Symbol==gene,]
   columns <- eval(parse(text=col))
   #print(columns)
   my.df <- df.child[,columns]
   
    if (length(columns)==1)
       {
       df.score <- my.df
       }
       # Use average of scores if more than one score is given
       else {df.score <- rowMeans(my.df)}
   
   my.seg <- protile.segmentation(gene,df.score,df.child$AA,half.size = as.numeric(size), 
                                  th.outlier = as.numeric(th1), th.tguh = as.numeric(th2))
     
   df.seg <- my.seg$segments     
   df.seg.filtered <- df.seg[df.seg$is.HS.site==TRUE,]
   df.seg.filtered$lenght <- df.seg.filtered$AA.end - df.seg.filtered$AA.start
   df.seg.filtered$Gene <- rep(gene,dim(df.seg.filtered)[1])
   
   df.score.aa <- data.frame('aver.score'=df.score,'aa.pos'=df.child$AA,'gene.symbol'=df.child$Symbol)
   v.filtered <- my.seg$v.filtered
   is.efficient <- my.seg$is.efficient
   df.score.aa <- cbind(df.score.aa[order(df.score.aa$aa.pos),],is.efficient,v.filtered)
    
   df.est.aa <- data.frame('aa.est'=my.seg$AA.est,'gene.symbol'=rep(gene,length(my.seg$AA.est)))
   
   seg.df = df.seg.filtered[,c(-1,-2)]
   score.df=df.score.aa 
   est.df=df.est.aa
   # Save table files in csv
   write.csv(seg.df,file=paste(gene,'_Segments.csv',sep=''),row.names=FALSE)
   write.csv(score.df,file=paste(gene,'_Score.csv',sep=''),row.names=FALSE)
   write.csv(est.df,file=paste(gene,'_Est.csv',sep=''),row.names=FALSE) 
   #return(list(seg.df = df.seg.filtered[,c(-1,-2)],score.df=df.score.aa,est.df=df.est.aa))
}

## Call HS region with user input parameters
CallHSRegion(args[1],args[2],args[3],args[4],args[5],args[6])
      