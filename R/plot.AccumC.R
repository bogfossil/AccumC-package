plot.AccumC <-
function(x,                 # an object of class 'AccumC'   
                        plot=T,            # make a simple plot in the current graphics window
                        bigplot=F,         # make a pdf showing all the calibrated dates, the age PDFs and the CAR
                        box_col='gray',    # the color of the bars on the CAR plot
                        prettytitle=NA,    # The title of the "bigplot"
                        filename = "AccumC_plot",     # directory location for "bigplot," defaults to main
                        ratelabel=T,       # add text labels indicating the rate for each bin
                        log_y_axis=F       # make the rate axis log
){
  
  if(class(x) != "AccumC") {stop('x does not have class "AccumC"', call. = F, domain = NA)} else {}
  
##--------------------------------------
## START Pull data from "AccumC" object
##--------------------------------------
  
  meta <- x$metadata
  bin_size <- x$metadata$bin_size
  calibration <- x$calibration
  rates <- x$rates
  prob_byone <- x$total_prob
  age_test <- x$age_test
  rate_test <- x$rate_test
  crate_test <- x$crate_test
  allpdfs <- x$all_pdfs
  AGES <- x$age_quantiles
  CRATES <- x$rate_quantiles
  binned <- x$binned
  maxLikeAge <- x$calibration$maxLikeAge
  
  
      DBS <- meta$bin_size
      DPT <- meta$depth
      AGE <- meta$age
      TOP <- round(min(calibration$maxLikeAge))
      ndbin <- round(DPT/DBS)
      corename <- meta$corename
#      if(is.na(corename) == T & bigplot == T){
#        while(file.exists(paste(corename))==F){
#          corename <- readline("What is the (case sensitive) folder name to save the bigplot?  \nIf you give the wrong answer, I'll ask you again.  ")
#        }
#          }else{}
      AGEM <- rates$age
      depths <- calibration$depths
      depthbins <- 0:(diff(range(depths))/bin_size)*bin_size + min(depths)
      age_bw <- meta$age_bw
       
      
##-------------------------------------
## END Pull data from "AccumC" object
##-------------------------------------
  

  

##-------------------------------------
## Start the little plot
##-------------------------------------

if(plot==T){
  
  par(mfrow=c(2,1), oma=c(4,4,4,4), mar=c(0,0,0,0), mgp=c(2, 1, 0))
  
  ### Age Model Plot
  plot(0.001*AGES[3,],DBS*0:(ndbin)+min(depths),
       xlim=c(0,AGE)*0.001, ylim=c(DPT,0), axes=F,
       xlab="", ylab="", font.lab=2, t="l", xpd=NA)
  
  axis(4)
  mtext(l=2.5, s=4, "Depth, cm", font=2)
  
  segments(0.001*AGES[1,], DBS*0:(ndbin)+min(depths), 
           0.001*AGES[5,], DBS*0:(ndbin)+min(depths),
           lty=3, lwd=1, xpd=NA)
  
  segments(0.001*AGES[2,], DBS*0:(ndbin)+min(depths), 
           0.001*AGES[4,], DBS*0:(ndbin)+min(depths),
           lty=1, lwd=4, col=box_col, xpd=NA)
  
  points(0.001*AGES[3,], DBS*0:(ndbin)+min(depths), pch=15, xpd=NA)
  
  if(ratelabel == T){
  text(0.001*AGES[5,-dim(AGES)[2]], DBS*1:ndbin-(DBS*0.5)+0.5*min(depths), 
       round(rates[,4], 3), cex=0.8, pos = 4)} else {}
  
  
  ### The carbon accumulation rate plot
  
  plot(AGEM*0.001, CRATES[3,], t="n",
       xlim=c(0,AGE)*0.001, ylim=range(CRATES[2:4,]), 
       axes=F, xlab="", ylab="", font.lab=2,
       log = ifelse(log_y_axis==T,"y",""))
  axis(2)
  mtext(l=2.5, s=2, "gC/m^2/a", font=2)
  axis(1)
  mtext(l=2.5, s=1, "Age, cal. ka", font=2)
  segments(0.001*AGEM, CRATES[1,], 0.001*AGEM, CRATES[5,], lwd=1, lty=3, xpd=NA)
  rect(0.001*AGES[3,-dim(AGES)[2]], c(CRATES[2,]), 0.001*AGES[3,-1], CRATES[4,], 
       col=box_col, xpd = NA)
  # points(0.001*AGEM,CRATES[3,], pch=15)
  segments(0.001*AGES[3,-dim(AGES)[2]], c(CRATES[3,]), 0.001*AGES[3,-1], CRATES[3,], lwd=2)
  
  
  
} else {}


##-------------------------------------
## End the little plot
##-------------------------------------

##-------------------------------------
## Start the BIG plot
##-------------------------------------

if(ndbin >= 100) {warning("Too many depth bins for bigplot. I set bigplot = FALSE", call. = F);
  bigplot = F} else {}
      
      
if(bigplot==T){
  
  
  ### Open a new pdf document  
  
  pdf(paste(filename,'.pdf',sep=""),w=6.5,h=9); par(oma=c(2,2,6,2))
  
  ### Set the layout for the page of plots
  
  layout(
    mat=array(c(1:(ncol(binned)-1)*0+ncol(binned),ncol(binned)+4,ncol(binned)+4,ncol(binned)+4,1:(ncol(binned)-1),ncol(binned)+1,ncol(binned)+2,ncol(binned)+3),
              dim=c(ncol(binned)+2,2)
    ),
    widths=c(1.25,3),
    heights=c(rep((1/(ncol(binned)-1)),(ncol(binned)-1)),0.3, 0.6, 0.6)
  )
  
  
  ### Plot PDFs for each depth bin
  
  for(i in 2:(ncol(binned))){
    par(mar=c(0,5,0,4)); plot(binned[,1],binned[,i],t='l',axes=F,ylab='',xlim=c(0,AGE)); box()
    mtext(s=2,l=1,paste(i-1),las=1)
  }
  
  
  title(outer=T,line=2,cex.main=2.5,
        main=ifelse(is.na(prettytitle),
                    paste(substr(toupper(corename),start=1,stop=1),substring(corename,first=2),sep=''),
                    prettytitle)
  )
  
  ### Plot all the calibrated dates  
  
  par(mar=c(0,4,0,0))
  plot(maxLikeAge*0.001,depths,ylim=rev(range(depthbins)),yaxs='i',
       ylab='Depth, cm',xlab='Age, cal. ka',font.lab=2,
       xpd=NA,cex.lab=1.2)
  abline(h=depthbins,lty=3)
  legend("topright", legend = "A", text.font = 2, cex = 2, bty = "n")
  
  #  par(mar=c(3.5,5,3,1));
  
  
  ### Lower Portion
  
  par(mar=c(0,5,3,4))
  
  ### combined probability plot
  
  plot(prob_byone$x*0.001, prob_byone$y, xlim=c(0, AGE)*0.001, axes=F, t="l", xlab="", ylab="")
  axis(2)
  mtext(l=2.5, s=2, "Probability", font=2)
  axis(3,at=0:30*2, line = 0.5)
  legend("topright", legend = "C", text.font = 2, cex = 2, bty = "n")
  
  ### Age Model Plot
  plot(0.001*AGES[3,],DBS*0:(ndbin)+min(depths),
       xlim=c(0,AGE)*0.001, ylim=c(DPT,0), axes=F,
       xlab="", ylab="", font.lab=2, t="l", xpd=NA)
  
  axis(4)
  mtext(l=2.5, s=4, "Depth, cm", font=2)
  
  segments(0.001*AGES[1,], DBS*0:(ndbin)+min(depths), 
           0.001*AGES[5,], DBS*0:(ndbin)+min(depths),
           lty=3, lwd=1, xpd=NA)
  
  segments(0.001*AGES[2,], DBS*0:(ndbin)+min(depths), 
           0.001*AGES[4,], DBS*0:(ndbin)+min(depths),
           lty=1, lwd=4, col=box_col, xpd=NA)
  
  points(0.001*AGES[3,], DBS*0:(ndbin)+min(depths), pch=15, xpd=NA)
  
  text(0.001*AGES[5,-dim(AGES)[2]], DBS*1:ndbin-(DBS*0.5)+0.5*min(depths), 
       round(rates[,4], 3), cex=0.8, pos = 4)
  
  legend("topright", legend = "D", text.font = 2, cex = 2, bty = "n")
  
  ### The carbon accumulation rate plot
  
  par(mar=c(4,5,0,4))
  plot(AGEM*0.001, CRATES[3,], t="n",
       xlim=c(0,AGE)*0.001, ylim=range(CRATES[2:4,]), 
       axes=F, xlab="", ylab="", font.lab=2,
       log = ifelse(log_y_axis==T,"y",""))
  axis(2)
  mtext(l=2.5, s=2, "gC/m^2/a", font=2)
  axis(1,at=0:30*2)
  mtext(l=2.5, s=1, "Age, cal. ka", font=2)
  segments(0.001*AGEM, CRATES[1,], 0.001*AGEM, CRATES[5,], lwd=1, lty=3, xpd=NA)
  rect(0.001*AGES[3,-dim(AGES)[2]], c(CRATES[2,]), 0.001*AGES[3,-1], CRATES[4,], 
       col=box_col, xpd = NA)
  # points(0.001*AGEM,CRATES[3,], pch=15)
  segments(0.001*AGES[3,-dim(AGES)[2]], c(CRATES[3,]), 0.001*AGES[3,-1], CRATES[3,], lwd=2)
  legend("topright", legend = "E", text.font = 2, cex = 2, bty = "n")
  
  
  
  ### Metadata box  
  
  
  par(mar=c(4.5,4,6,0)); 
  plot(0:1, 0:1 ,axes=F,xlab='',ylab='', pch=""); box()
  
  text(rep(1, 6), 6:1/7, round(c(DPT, min(AGES[3,]), max(AGES[3,]), meta$ndates, DBS, age_bw)),
       adj=1, font=2)
  text(rep(0, 6), 6:1/7, paste(c("Thickness (cm)", "Top Age", "Bottom Age", "n Dates", "Bin Size (cm)", "Age BW (yr)"),
                               ': ',sep=''), adj=0, font=3)
  
  legend("topright", legend = "F", text.font = 2, cex = 2, bty = "n")
  
  mtext("B", side = 3, line = -2, outer = T, adj = 0.9, cex = 1.5, font = 2)
  
  
  
  dev.off()
  
} else {}

##-------------------------------------
## End the BIG plot
##-------------------------------------

      
      
}
