## Source file for qPCR analysis

library(RDML)
library(qpcR)
library(dplyr)
library(scales)
library(ggplot2)
library(reshape2)


### Function definitions

## Read RDML file, fit models and extract relevant data

process.run <- function(input, infofile = NULL, opt = FALSE, exclude = NULL, fixbug = F) {

  rdml <- RDML$new(input)
  info <- rdml$AsTable()
  
  
  #rdml$AsDendrogram()
  amp <- rdml$GetFData(filter(info, sample.type == "unkn"))
  mods <- modlist(amp, 1
                  , model = l5
                  , check = "uni2"
                  #, checkpar = parKOD()
                  #, remove = "KOD"
                  , exclude = exclude
                  , labels = NULL
                  , norm = TRUE
                  , baseline = "mean"
                  , basecyc = 1:8
                  , basefac = 1
                  , smooth = "savgol"
                  #, smoothPAR = NULL
                  , opt = opt
                  , optPAR = list(crit = "weights")
                  , verbose = TRUE
  )
  
  
  fits <- pcrbatch(mods, cyc = 1, fluo = NULL, 
                   methods = c("sigfit","cm3"),
                   type = "Cy0",
                   #labels = NULL,
                   #group = NULL,
                   #names = c("group", "first"),
                   plot = TRUE, 
                   verbose = TRUE)
  
  # workaround for bug in qpcR
  if(fixbug){
    for(i in 1:ncol(fits)) {
       if(!is.null(dim(fits[,i]))) {
         fits[,i] <- fits[,i][,1]
       }
    }
  }


  # make data frame
  
  fits.df <- as.data.frame(fits[,-1], stringsAsFactors=F)
  fits.df.t <- data.table::transpose(fits.df)
  colnames(fits.df.t) <- fits[,1]
  fits.df.t$fdata.name <- colnames(fits[,-1])
  
  info.df <- as.data.frame(rdml$AsTable(), stringsAsFactors=F)
    full <- inner_join(fits.df.t, info.df)
  
  d0 <- data.frame(full$sample,
                   full$target,
                   as.numeric(full$cm3.D0),
                   as.numeric(full$sig.Cy0),
                   as.numeric(full$sig.eff),
                   stringsAsFactors=F)
  
  colnames(d0) = c("sample", "target", "d0", "cy0", "eff")
  
  if(!is.null(infofile)) {
    sample.info <- read.table(infofile, sep="\t", header = T)
    group <- sample.info[match(d0$sample, sample.info[,1]), 4, drop=F]
    d0$group <- group[,1]
  }
  
  plot(mods, type="single", add=F, col = as.factor(sapply(mods, function(x) strsplit(x$names, "_")[[1]][5])))
  return(d0)
  
}


# function for geometric mean
# stolen from https://stackoverflow.com/questions/2602583

gm_mean <- function(x, na.rm=TRUE) {
  
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  
}

## normalize expression.
##     one reference method: d0(goi)/d0(ref)
##     multiple references: d0(goi)/gm_mean(d0(refs))

normalize <- function(input, ref) {
  
  # TODO: handle technical replicates
  
  if(length(ref) == 1) {
    normalized <- input %>% group_by(sample) %>% mutate(d0.norm = d0/d0[ target==ref ]) # gotta love dplyr
  }
  
  if(length(ref) > 1) {
    normalized <- input %>% group_by(sample) %>% mutate(norm.f = gm_mean(d0[ target %in% ref ]),
                                                        d0.norm = d0/norm.f)
  }
  
  return(normalized)
}


## normalize using Cy0 instead of D0

normalize.cy0 <- function(input, ref) {
  
  # TODO: handle technical replicates
  
  if(length(ref) == 1) {
    normalized <- input %>% group_by(sample) %>% mutate(cy0.norm = cy0[ target==ref ] - cy0) 
  }
  
  if(length(ref) > 1) {
    normalized <- input %>% group_by(sample) %>% mutate(norm.f = gm_mean(cy0[ target %in% ref ]),
                                                        cy0.norm = norm.f - cy0)
  }
  
  return(normalized)
}


## summarize expression by group; column name for grouping variable is hardcoded -
##   easier than messing around with dplyr's NSE

summarize.ex <- function(input) {
  
  summa <- input %>%
    group_by(group, target) %>%
    summarize(mean.ex = mean(d0.norm),
              sd.ex = sd(d0.norm),
              sem.ex = sd.ex/sqrt(n()),
              n = n())
  
  return(summa)
}


## summarize expression using Cy0

summarize.ex.cy0 <- function(input) {
  
  summa <- input %>%
    group_by(group, target) %>%
    summarize(mean.ex = mean(cy0.norm),
              sd.ex = sd(cy0.norm),
              sem.ex = sd.ex/sqrt(n()),
              n = n())
  
  return(summa)
}

### Function definitions end

message("Loaded functions: process.run, gm_mean, normalize, normalize.cy0, summarize.ex, summarize.ex.cy0")
