find.Xzone.Joe <- function(output, 
         BS.dive=F, plot=TRUE, plot.legend=TRUE, cex.numbers=1.3, cex.axis=1.3, cex.lab=1.3, lwd=2, dev.new=TRUE, plot.xzone=TRUE, plot.truth=TRUE, plot.bsm=TRUE, plot.title=TRUE, draw.numbers=TRUE){
  
  ## ARGUMENTS TO FUNCTION - use args(find.Xzone) to see the defaults

  # output: the object holding the output from get.BSMdive 
  # BS.dive: is a logical argument indicating whether the dive data have already been abstracted by the broken-stick model
  # plot: is a logical argument and indicates whether you want your dive to be plotted
  # plot.legend: is a logical argument and indicates whether you want ta legend to be plotted
  # cex.numbers: is a real number and indicates the size of the numbers indicating the BSM points
  # cex.axis: is a real number and indicates the size of the numbers in the axes
  # cex.lab: is a real number and indicates the size of the axis labels 
  # lwd: is a real number and indicates the width of the line of the dive zone  
  # dev.new: is a logical argument and indicates whether you want your dive to be plotted in a new window 
  # plot.xzone: is a logical argument and indicates whether you want the dive zone to be plotted
  # plot.truth: is a logical argument and indicates whether you want the detailed dive to be plotted. this only applies when you have detailed dive data. 
  # plot.bsm: is a logical argument and indicates whether you want the broken-stick abstracted dive to be plotted
  # plot.title: is a logical argument and indicates whether you want a title to be plotted
  # draw.numbers: is a logical argument and indicates whether you want numbers indicating the order of the points selected by the broken-stick be plotted
  
  depth.bin.table <- output$depth.bin.table
  ds <- output$ds
  BS.time.depth <- output$BS.time.depth  
  s <- output$s
  breakpoints <- output$breakpoints 
  BS.time.depth <- output$BS.time.depth  
  time.depth <- as.data.frame(output$time.depth) 
  BS.order <- output$BS.order 
  sampletime.steps <- output$sampletime.steps 
  sampletime.res <- output$sampletime.res
  divetable <- output$divetable 
  dn <- output$dn 
  step.size.vec <- output$step.size.vec 
  R.vec <- output$R.vec # objects from both detailed and abstracted dives have this, but note R5 from bsmtest is huge, whereas from dtest its precise
  R5.trans <- output$R5.trans 
  sample.index <- output$sample.index 
  dloc <- output$dloc 
  seg.length <- output$seg.length 
  step.size.vec <- output$step.size.vec
  
  ##############################################################
  # Code up the dive zone...
  ##############################################################
  
  # ii is dimension y (of [x,y,z]) of the ds array, normally y=5 or more
  # NB. for 4 internal points, I will have 6 sets of coordinates and 5 line segments
  ii <- ds.col <- ncol(ds)
  
  # tt is the sequence of times through the dive where you want to *generate* upper and lower bounds
  # it has to start with 1 not 0 because its also an index. 
  
  tt <- seq(1, nrow(BS.time.depth)) 
  depths <- BS.time.depth$depths
  
  R.index <- Resid.index <- ii-1 
  
  # array of dim(times, ii, 2) to store upper and lower locations for each segment
  tempX.lines <- array(0, dim=c(length(tt),ii,2)); #dim(tempX.lines)
  allX.lines <- array(0, dim=c(length(tt),ii,2)); #dim(allX.lines)
  
  # array to store final upper and lower limits (2 columns) of exclusion zone for each time point, length(tt) rows
  Xzones <- array(0, dim=c(length(tt),2*ii)); #dim(Xzones)
  
  dimnames(Xzones) <- list(NULL,paste(c("tminR", "tmaxR"), rep(1:ds.col, each=2), sep="")) # tminR is the upper limit, tmaxR is the lower limit (counter-intuitive, I know)
  
  calc.Xzones <- function(){
    
    # create two vectors, UPP and LOW, of length ii, to hold the coordinates of the upper and lower plausible bounds for each segment 
    UPP <- rep(0, ii) # deep
    LOW <- rep(0, ii) # shallow
    BSM <- matrix(data=NA, nrow=length(tt), ncol=ii)
    # create a vector, m, of length ii, to hold the slope of each segment considered
    m <- rep(0, ii)
    r <- rep(0, ii)
    
    for (t in tt){
      # run through each sample time point (tt)
      
      for (j in 1:(ii)){ # run through each iteration/segment of the model
        # j will run through all segments for each time step t
        
        # start with i being 1, ie the first element of v
        i <- 1
        
        # let v be the jth row of the first "table" in the 3D s array. this will be a vector 5 elements long
        v <- s[j,,1]
        
        # populate the ii+1th element with the number of samples
        v[ii+1] <- length(tt) 
        
        # take the 1st element from the jth row of the s array and sort it in ascending order 
        v <- sort(v)
        
        # as long as the ith element of v is smaller than the chosen t then move to the next element until you find one that is bigger than t
        while (v[i] < t) {i <- i+1}
        
        if (v[i] == t) { # AT BREAKPOINTS
          
          bin.row.index <- which(depths[t] >= depth.bin.table[,1] & depths[t] < depth.bin.table[,2])
          step.size <- depth.bin.table[bin.row.index, 3] 
          
          # I only need to add (not subtract as well) the step size to the depth because the depth that is returned, depths[t], is the bottom (shallow end) of the bin.
          UPP[j] <- min(depths[t] + step.size, max(depths) + step.size.vec[sort(BS.order)==1])  # deeper than profile
          LOW[j] <- depths[t]               # shallower than profile
          if(t %in% v[1:j]){ BSM[t,j] <- depths[t] } 
          # this still leaves some NAs in the BSM matrix, for some reason...
          
        } else { # this bracket closes the if(v[i] == t) statement
          
          # the ith element of v will be the first element that is bigger than t, this will be the upper x bound of the time bin that t belongs to
          up <- v[i]
          
          # therefore, the the (i-1)th element of v will be the element that is immediately smaller than t, 
          # this will be the lower x bound of the time bin that t belongs to
          lo <- v[i-1]
          
          # m gives the slope of the segment in question
          m[j] <- (depths[up]-depths[lo])/(up-lo)
          
          # r gives the vertical distance the exclusion zone should be drawn at, above and below the point in the profile, that is being considered
          r[j] <- max(R.vec[j], min(step.size.vec)) 
          
          seg.b <- depths[lo] - m[j]*lo
          
          # this loop gives me depths along the line segment and the equation of that line so that I can project h above and below that line. the reason I need to do this is that at the moment my simulated dive is on a much coarser timescale that the time points I am "sampling" so that I can draw the exclusion lines
          
          BSM[t,j] <- m[j]*t + seg.b
          
          # This is meant to define the actual values of the upper and lower limits to draw on the graph 
          UPP[j] <- BSM[t,j] + r[j]  # deeper than profile
          LOW[j] <- BSM[t,j] - r[j]  # shallower than profile
          
        } # this bracket closes the if(v[i] == t) else statement
        
      } # this bracket closes the for (j in 1:ii) loop - segments/iterations
      
      # zone limits before maximum depth and surface limits are in place		
      allX.lines[which(tt==t),,1] <- LOW # shallower
      allX.lines[which(tt==t),,2] <- UPP # deeper
      
      # work out the minimum width of the dive zone along the broken-stick profile
      interp.up <- interp.lo <- 0
      s.dloc <- sort(dloc)
      for (i in 1:(length(dloc)-1)){	
        from.lo <- BS.time.depth[s.dloc[i],2]
        to.lo <- BS.time.depth[s.dloc[i+1],2]	
        from.up <- BS.time.depth[s.dloc[i],2]+step.size.vec[i]
        to.up <- BS.time.depth[s.dloc[i+1],2]+step.size.vec[i+1]	
        if (BS.dive==T) {length.out <- seg.length[i]+1} else {length.out <- seg.length[i]}
        interp.lo[s.dloc[i]:s.dloc[i+1]] <- seq(from=from.lo, to=to.lo, length.out=length.out)
        interp.up[s.dloc[i]:s.dloc[i+1]] <- seq(from=from.up, to=to.up, length.out=length.out)
      }
      
      # loop through each iteration and implement surface, maximum depth and other limits on current dive zone limits	
      for (q in 1:ii){
        
        # zone limits including surface and maximum depth limits
        Xzones[t,(2*q)] <- tempX.lines[t,q,1] <- max(0, max(allX.lines[t,c(1:(ii-1)),1])) # LOW:  because depths are positive, the max is the deepest value of the shallow boundary	
        min.X <- min(allX.lines[t,c(1:(ii-1)),2])
        Xzones[t,(2*q-1)] <- tempX.lines[t,q,2] <- max(min(max(depths) + step.size.vec[sort(BS.order)==1], min.X), interp.up[t]) # UPP: because depths are positive, the min is the shallowest value of the deep bounday. also, you want the zone to be at least the bin width of the break point. this avoids zero width at breakpoints
        
      } # q loop
      
    } # this bracket closes the for (t in tt) loop
    
    Xout <- list(Xzones, BSM, allX.lines)
    
  } # close function calc.Xzones() here and check if any points in the true-astracted profile fall outside the Xzone. if yes, rerun the function with Resid.index <- ii-1	
  
  Xout <- calc.Xzones();
  Xzones <- Xout[[1]]; BSM <- Xout[[2]]; allX.lines <- Xout[[3]]
  
  ####### adding depths bins where breakpoints occur
  
  # Incorporating Rmax in dive zone 
  t.Rmax <- s[R.index+1,R.index+1,1]
  Rmax.LO <- Xzones[t.Rmax, (2*(ii-1)-1)]==BS.time.depth[s[R.index+1,R.index+1,1],2] # find the row in the Xzone matrix where the lower boundary touches the dive path
  Rmax.UP <- Xzones[t.Rmax, (2*(ii-1))]==BS.time.depth[s[R.index+1,R.index+1,1],2] # find the row in the Xzone matrix where the upper boundary touches the dive path
  
  if(Rmax.LO == TRUE){ # if its not 0, ie TRUE outcome to trial, it means this is the side of the zone the last BS point touches, so dont change anything in the LOWER BOUNDARY and bring the UPPER BOUNDARY up to the depth + step.size for its relevant depth bin
    
    Rmax.UP <- BS.time.depth[t.Rmax,2] + step.size.vec[sort(BS.order)==R.index]
    Xzones[t.Rmax, (2*(ii-1)-1)] <- Rmax.UP
    
  } else { # if its not 0 it means this is the side of the zone the last BS point touches, so bring the LOWER BOUNDARY down to the shallow end of its relevant depth bin, and increase the UPPER BOUNDARY to the deep end of its relevant depth bin
    
    Rmax.LO <- BS.time.depth[s[R.index+1,R.index+1,1],2]
    Rmax.UP <- BS.time.depth[s[R.index+1,R.index+1,1],2] + step.size.vec[sort(BS.order)==R.index]
    
    # sometimes you will have two points in a dive that touch a boundary. to stop this causing problems I'm going to make it that one of those points gets picked randomly. it only really serves to create the boundary so it should be ok to do it like this. 
    
    if(length(Rmax.LO)>1 | length(Rmax.UP)>1){
      Rmax.LO <- Rmax.LO[1]
      Rmax.UP <- Rmax.UP[1]
    }
    
    Xzones[t.Rmax, (2*(ii-1))] <- Rmax.LO
    Xzones[t.Rmax, (2*(ii-1)-1)] <- Rmax.UP
    
  } # closes if(Rmax.LO == TRUE)
  
  ###### work out the dive zone index
  maxdep <- max(BS.time.depth$depths)
  sampletime.res <- output$sampletime.res # should be 4, 8, 16...
  tmax <- max(BS.time.depth$times)*60/sampletime.res 
  dz.height <- vector("numeric", length=nrow(BS.time.depth))
  
  for (i in 1:nrow(BS.time.depth)){
    dz.height[i] <- Xzones[i,(2*R.index)-1]-Xzones[i,2*R.index] # Upper (deep) limit of dive zone MINUS Lower (shallow) limit of dive zone	
  }
  ######
  
  dz.index <- sum(dz.height)/(maxdep*tmax)
  
  ########################################################
  # Plotting
  ########################################################
  bs.points <- data.frame(time=R.index+2, depth=2)
  
  for(i in 1:(R.index+1)){
        bs.points[i,2] <- BS.time.depth[s[R.index+1,i,1],2]
    bs.points[i,1] <- BS.time.depth[s[R.index+1,i,1],1]
  }
  
  bs.points[11,2] <- BS.time.depth[s[1,1,2],2]
  bs.points[11,1] <- BS.time.depth[s[1,1,2],1]
  
  bs.points <- bs.points[order(bs.points$time),]
  
  if (plot==TRUE){
    
    if (dev.new==TRUE){
      dev.new(); 
    }
    
    par(mar=c(5,5,1,1))
    
    if (plot.truth==TRUE){
      plot(time.depth$times, -time.depth$depths, ylim=c(-(max(time.depth$depths)+0.3*max(time.depth$depths)), 0.1*max(time.depth$depths)), ylab="Depth (m)", xlab="Time (min)", col=2, type="l", lty=3, lwd=lwd, cex.axis=cex.axis, cex.lab=cex.lab)
      
    } else {
      
      plot(BS.time.depth$times, -BS.time.depth$depths, ylim=c(-(max(BS.time.depth$depths)+0.3*max(BS.time.depth$depths)), 0.1*max(BS.time.depth$depths)), ylab="Depth (m)", xlab="Time (min)", col=2, type="n", lty=3, cex.axis=cex.axis, cex.lab=cex.lab, lwd=lwd)
      
    } # this bracket closes the if(plot.truth==TRUE) statement			
    if (plot.bsm==TRUE){
      bs.points <- matrix(nrow=R.index+2, ncol=2)

            for(i in 1:(R.index+1)){
        
        segments(x0=BS.time.depth[s[R.index+1,i,1],1], y0=-BS.time.depth[s[R.index+1,i,1],2], x1=BS.time.depth[s[R.index+1,i,2],1], y1=-BS.time.depth[s[R.index+1,i,2],2], lwd=lwd-1)
        bs.points[i,2] <- BS.time.depth[s[R.index+1,i,1],2]
        bs.points[i,1] <- BS.time.depth[s[R.index+1,i,1],1]
      }
      
      bs.points[11,2] <- BS.time.depth[s[1,1,2],2]
      bs.points[11,1] <- BS.time.depth[s[1,1,2],1]
      
      if (draw.numbers==TRUE){
        
        text(x=BS.time.depth[s[1,1,1],1],y=-(BS.time.depth[s[1,1,1],2]-0.1*max(BS.time.depth$depths)),paste("0"), cex=cex.numbers) 
        text(x=BS.time.depth[s[1,1,2],1],y=-(BS.time.depth[s[1,1,2],2]-0.1*max(BS.time.depth$depths)),paste(R.index+1), cex=cex.numbers) 
        
        for(i in 2:(R.index+1)){
          
          text(x=BS.time.depth[s[i,i,1],1],y=-(BS.time.depth[s[i,i,1],2]+0.1*max(BS.time.depth$depths)),paste(i-1), cex=cex.numbers) 
          
        }
        
      } # this bracket closes the if(draw.numbers==TRUE) statement
      
    } # this bracket closes the if(plot.bsm==TRUE) statement
    
    if(plot.xzone==TRUE){
      lines(BS.time.depth[c(1:length(tt)),1], -Xzones[,(R.index*2)-1], col="orange", lwd=lwd, type="l")
      lines(BS.time.depth[c(1:length(tt)),1], -Xzones[,R.index*2], col="orange", lwd=lwd, type="l")
    }
    
    if (plot.legend == TRUE){
      if (plot.xzone==TRUE & plot.truth==TRUE){
        legend("bottomleft", legend=c(paste("Dive zone ", "(DZI: ", round(dz.index, digits=2), ", R", Resid.index, ": ", round(R.vec[Resid.index]), "m", ")", " ", sep=""), "Detailed profile", "BSM profile"), bty="n", lty=c(1,2,1), col=c("orange", "red", 1), lwd=c(lwd,lwd-1,lwd-1), cex=cex.lab)
      } else {
        if (plot.truth==FALSE){
          legend("bottomleft", legend=c(paste("Dive zone ", "(DZI: ", round(dz.index, digits=2), ", R", Resid.index, ": ", round(R.vec[Resid.index]), "m", ")", " ", sep=""), "BSM profile"), bty="n", lty=c(1,1), col=c("orange", 1), lwd=c(lwd,lwd-1), cex=cex.lab)
        } else {
          if (plot.xzone==FALSE){
            legend("bottomleft", legend=c("Detailed profile", paste("BSM profile ", "(DZI: ", round(dz.index, digits=2), ", R", Resid.index, ": ", round(R.vec[Resid.index]), "m", ")", " ", sep="")), bty="n", lty=c(2,1), col=c("red", 1), lwd=c(lwd-1,lwd-1), cex=cex.lab)
          } # closes if
        } # closes else
      }
    } # plot.legend bracket
    
    if(plot.title == TRUE){
      if(BS.dive==T){
        title(main=paste(paste("Dive",dn,sep=" "),paste("(",divetable$ref[dn],")",sep=""),sep=" "))		
      } else {
        title(main=paste(paste("Dive",dn,sep=" "),paste("(",divetable$ref[dn],")",sep=""),sep=" "))	
      }
      
    } # closes if(plot.title==TRUE)
    
  }
  
  out <- list(bs.points=bs.points, tt=tt, allX.lines=allX.lines, BS.time.depth=BS.time.depth, Xzones=Xzones, R.index=R.index, Resid.index=Resid.index, R.vec=R.vec, dive.number=dn, DZI=dz.index) 
  return(out)
  
}
