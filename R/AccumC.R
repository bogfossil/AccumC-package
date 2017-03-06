AccumC <-
function(x = imv,           # a data frame where name(x) == c("labcode", "depth", "bp", "std", "afbd")
                   corename=NA,       # the name of the directory and filenames for "corename.dat" and "corename.afbd"
                   bin_size=NA,       # the tickness of the depth bins in centimeters,
                   age_bw = NA,       # the size of the age bandwith in years
                   kernel = "r",      # Kernel shape for PDF smoothing, as in stats::density()
                   calcurve=intcal13, # can be intcal13 or shcal13 or something else you have
                   it_power=17,       # log[2] of the number of iterations for Monte Carlo simulation
                   pct_C=0.48         # fraction of carbon per unit organic matter (translates afbd into carbon)
                   ){
  
data(intcal13)
data(shcal13)

it <- 2^it_power

##-------------------------------------
## Start import data
##-------------------------------------
 if(sum(c("labcode", "depth", "bp", "std", "afbd") %in% names(imv))==5){
   datfile <- x
   }else{ 
     
     while(file.exists(paste(corename,"/",corename,".dat",sep=""))==F){
       corename <- readline("What is the (case sensitive) core name?  \nIf you give the wrong answer, I'll ask you again.  ")
     }
     
     datfile <- read.table(paste(corename,"/",corename,".dat",sep=""),header=T,sep=",")
   }

datfile <- datfile[order(datfile$bp),]
labcodes <- datfile$labcode

##-------------------------------------
## End import data
##-------------------------------------

##-------------------------------------
## Start map dates to cal curve
##-------------------------------------

allpdfs <- list(NA)
for(i in 1:length(labcodes)) {
  bp <- (-5*datfile$std[i]+datfile$bp[i]):(5*datfile$std[i]+datfile$bp[i]) 
  p <- dnorm(bp,datfile$bp[i],datfile$std[i]) 
  date <- data.frame(bp,p)
  
  cal <- merge(calcurve,date,by='bp')
  cal <- data.frame(cal$cal,cal$p)
  cal[,2] <- cal[,2]/sum(cal[,2])
  allpdfs[[i]] <- cal[order(cal[,1]),]
}

for(i in 1:length(allpdfs)){
  names(allpdfs[[i]]) <- c("age",labcodes[i])
}


merged <- allpdfs[[1]]  # Put all the dates onto a single age scale
for(i in 2:length(allpdfs)){
  merged <- merge(merged,allpdfs[[i]],all=T)
}
merged[is.na(merged)] <- 0 #make NAs in the PDF array into zeroes



depths <- datfile$depth #These are the depths of each date

if(file.exists(paste(corename,"/",corename,".afbd",sep=""))==T) {
  bulkDen <- read.table(paste(corename,"/",corename,".afbd",sep=""),header=T,sep=",")
  datfile$afbd <- approx(bulkDen,xout=depths,rule=2)$y
} else {}

maxLikeAge <- NA
for(i in 2:dim(merged)[2]) {
  maxLikeAge[i-1] <- merged[,1][merged[,i]==max(merged[,i])][1] 
}

##-------------------------------------
## End map dates to cal curve
##-------------------------------------

##-------------------------------------
## Start merging age PDFs
##-------------------------------------

# Bin the dates by depth -- 3 (or so) dates per depth bin 
# (make it an integer, always round down, hence "floor")

  if(is.na(bin_size)==T) {
    bin_size <- diff(range(depths))/floor(length(labcodes)/3)
  } else {
    nbins <- round(diff(range(depths))/as.integer(bin_size))
    bin_size <- diff(range(depths))/nbins
  }

depthbins <- 0:(diff(range(depths))/bin_size)*bin_size + min(depths)
binned <- array(dim=c(nrow(merged),length(depthbins))) #make an empty array with (n-years) rows and (n-depthbins) columns
binned[,1] <- merged[,1] #make the first column years
merged <- merged[,-1] #remove the "years" column from the PDF array

for(i in 2:length(depthbins)){  #make a mean PDF for each bin
  
  xx <- merged[,depths>=depthbins[i-1]&depths<=depthbins[i]] #assign the bin to a variable
  
  if(mode(xx)=="numeric") {binned[,i] <- xx} else {binned[,i] <- apply(xx,1,mean)} 
  #if there's only one date in a bin, it will be "numeric" else it's a "list"
  
}
binned[is.na(binned)] <- 0
binned[is.nan(binned)] <- 0

### Combine all the probabilities

prob <- apply(binned[,-1],1,sum)/(ncol(binned)-1)
prob[is.na(prob)] <- 0
prob <- data.frame(binned[,1], prob)

##-------------------------------------
## End merging age PDFs
##-------------------------------------


# Tabluate Metadata
calibration <- data.frame(labcodes,depths,maxLikeAge)
meta <- data.frame(bin_size,
                   length(labcodes), 
                   round(max(binned[,1])), 
                   round(bin_size*(ncol(binned)-1)),
                   corename)
names(meta) <- c("bin_size","ndates", "age", "depth", "corename")

# Carbon Data
afbd <- data.frame(maxLikeAge,datfile$afbd) 
afbd <- afbd[order(maxLikeAge),]
names(afbd) <- c("depth", "afbd")


# -----------------------------------------------------------------------
# Calculate the Age Model and C-accumulation rates WITH UNCERTAINTIES!!
# 
# Start 
# -----------------------------------------------------------------------

DBS <- meta$bin_size

DPT <- meta$depth
AGE <- meta$age
TOP <- round(min(calibration$maxLikeAge))

ndbin <- round((DPT/DBS))

if(is.na(age_bw)==T) {
  age_bw <- diff(c(TOP,AGE))/length(labcodes)
} else {
  age_bw <- as.integer(age_bw)
}

prob_byone <- density(sample(x = prob[,1], size = it, replace = T, prob = prob[,2]),
                      cut = 0, bw = age_bw, kernel = kernel)
prob_byone <- approx(prob_byone, xout=TOP:AGE)
prob_byone$y[is.na(prob_byone$y)] <- 0

age_test <- array(
  sample(TOP:AGE, replace=T, size=(ndbin+1)*it, prob=prob_byone$y),
  dim=c(ndbin+1, it)
)

max_age <- apply(age_test, 2, min)-apply(age_test, 2, max)

age_test <- apply(age_test, 2, sort)
age_test <- array(age_test, dim=c(ndbin+1, ncol(age_test)/3, 3) )
age_test <- apply(age_test, 1:2, median)

#age_test[1, ] <- sample(allpdfs[[1]][,1], size = dim(age_test)[2], replace = T, prob = allpdfs[[1]][,2])
#age_test[dim(age_test)[1], ] <- sample(allpdfs[[length(allpdfs)]][,1], size = dim(age_test)[2], replace = T, prob = allpdfs[[length(allpdfs)]][,2])
  
#age_test <- age_test[, age_test[1, ] <= max(allpdfs[[1]])]
#age_test <- age_test[, age_test[dim(age_test)[1], ] >= min(allpdfs[[length(allpdfs)]])]
 
rate_test <- DBS/apply(age_test, 2, diff)

##-------------------------------------
## Start get rid of unreasonable rates
##-------------------------------------

#avg_rate <- (ndbin*DBS)/apply(age_test, 2, max)  ## find the average rate for each iteration
#avg_rate_range <- (ndbin*DBS)/(range(allpdfs[calibration$maxLikeAge==max(calibration$maxLikeAge)][[1]]$age)-range(allpdfs[calibration$maxLikeAge==min(calibration$maxLikeAge)][[1]]$age))
       ## find out the range of bottom dates, therefore the range of average sed rates
#avg_top_range <- (ndbin*DBS)/(range(allpdfs[calibration$maxLikeAge==max(calibration$maxLikeAge)][[1]]$age)-range(allpdfs[calibration$maxLikeAge==min(calibration$maxLikeAge)][[1]]$age))
       ## find out the range of bottom dates, therefore the range of average sed rates

# rate_test <- rate_test[,avg_rate<=avg_rate_range[1]&avg_rate>=avg_rate_range[2]]    ## remove any tests that have a bottom date younger than the oldest calibrated date range
# rate_test <- rate_test[,avg_rate<=avg_rate_range[1]&avg_rate>=avg_rate_range[2]]    ## remove any tests that have a top date older than the youngest calibrated date range
rate_test <- rate_test[,apply(is.finite(rate_test), 2, sum)==nrow(rate_test)]       ## remove any tests that include infinite sed rates
# rate_test <- rate_test*(1/median(as.numeric(rate_test)))*(diff(range(calibration$depths))/diff(range(calibration$maxLikeAge)))

##-------------------------------------
## End get rid of unreasonable rates
##-------------------------------------


AFBD <- matrix(approx(afbd, n=ndbin)$y, byrow=F, nrow = nrow(rate_test), ncol = ncol(rate_test))
crate_test <-rate_test*AFBD*pct_C*10000
CRATES <- apply(crate_test, 1, quantile, na.rm=T, probs=c(0.125, 0.25, 0.5, 0.75, 0.875))
AGES <- apply(apply(age_test, 1:2, median), 1, quantile, na.rm=T, probs=c(0.125, 0.25, 0.5, 0.75, 0.875))
AGEM <- apply(age_test, 1, median)[-1] - 0.5*diff(apply(age_test, 1, median))

crate_byone <- data.frame(prob_byone$x, 
                          diff(range(depths))*prob_byone$y*approx(afbd, n=length(prob_byone$x))$y*10000)
names(crate_byone) <- c("age", "crate")

borders <- apply(apply(age_test, 1:2, median),1,median)
rates <- data.frame(AGEM,borders[-1], borders[-length(borders)],apply(rate_test, 1, median),CRATES[3,], apply(AFBD, 1, median))
names(rates) <- c("age","bottom_age","top_age","sed. rate (cm/a)","C accum. rate (g/m^2/a)", "afbd (g/cc)")

# -----------------------------------------------------------------------
# Calculate the Age Model and C-accumulation rates WITH UNCERTAINTIES!!
# 
# End 
# -----------------------------------------------------------------------


## Sew up all the output into a "list" object
meta$age_bw <- age_bw
output <- list(meta, calibration, rates, prob_byone, AGES[,-1]-t(apply(AGES, 1, diff)*0.5), age_test, rate_test, crate_test, allpdfs, AGES, CRATES, binned)
names(output) <- c("metadata","calibrations","rates", "total_prob", "age_mid_quantiles", "age_test", "rate_test", "crate_test", "all_pdfs", "age_quantiles", "rate_quantiles", "binned_pdfs")

class(output) <- "AccumC"

return(output)

}
