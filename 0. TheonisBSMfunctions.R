################################################################################################
# THIS FUNCTION WORKS ON ABSTRACTED (BROKEN-STICK) OR DETAILED TIME-DEPTH DIVE DATA FROM SEALS
# COLLECTED USING CTD-SRDLs OR SRDLs
# T. Photopoulou, P. Lovell and J. Matthiopoulos
# Last edits 20141201
################################################################################################

################ get.BSMdive
# This function does two things, depending on the kind of dive data you give it.
# a) If you give it a vector of depths sampled at regular intervals, it carries out abstraction using the broken-stick model
# b) If you give it a dive that has already been abstracted on-board an SRDL, it works out the order in which the points were added to the abstracted profile 

get.BSMdive <- function(BS.dive=FALSE, detailed.depth, res=4, divetable, breakpoints=4, dn=1, plot=TRUE, dev.new=TRUE, plot.truth=TRUE, plot.bsm=TRUE, cex.axis=1.3, cex.lab=1.3, draw.numbers=TRUE, draw.lines=TRUE){
  
  ## ARGUMENTS TO FUNCTION - use args(get.BSMdive) to see the defaults
  # BS.dive: is a logical argument indicating whether the dive data have already been abstracted by the broken-stick model
  # detailed.depth: is a vector of depths and this argument should be provided when you have high resolution data
  # res: takes an integer and is the resolution of the high resolution data in seconds, its default is 4 sec
  # divetable: is a dataframe and should be provided when you have already abstracted dive data. this should have column titles the same as the SMRU divetable's 
  # breakpoints: is an integer and is the number of breakpoints you want in your abstracted dive. the default is 4
  # dn: is an integer indicates the dive number, essentially the row index in the divetable for the dive you want to abstract
  # plot: is a logical argument and indicates whether you want your dive to be plotted
  # dev.new: is a logical argument and indicates whether you want your dive to be plotted in a new window
  # plot.truth: is a logical argument and indicates whether you want the detailed dive to be plotted. this only applies when you have detailed dive data.
  # plot.bsm: is a logical argument and indicates whether you want the broken-stick abstracted dive to be plotted
  # cex.axis: is a real number and indicates the size of the axis
  # cex.lab: is a real number and indicates the size of the axis labels
  # draw.numbers: is a logical argument and indicates whether you want numbers indicating the order of the points selected by the broken-stick be plotted
  # draw.lines: is a logical argument and indicates whether you want lines to be plotted, one horizontal line at max.depth and one vertical one at the time of max.depth
  
  # make sure your R session uses UTC for all your datetimes unless some other time zone is required
  Sys.setenv(TZ='UTC')
  
  # you want these bins to be open at the top [0, 1575). depths are divided by 10 because the tag deals in decimeters
  depth.bins <- c(0, 1575, 1612.5, 4762.5, 4837.5, 11137.5, 11287.5, 23887.5)/10
  depth.step <- c(25, 37.5, 50, 75, 100, 150, 200)/10
  
  bin.table.colnames <- list("shallow", "deep", "step")
  depth.bin.table <- matrix(data=0, nrow=7, ncol=3, dimnames=list(NULL,bin.table.colnames))
  depth.bin.table[,1] <- depth.bins[1:7]
  depth.bin.table[,2] <- depth.bins[2:8]
  depth.bin.table[,3] <- depth.step 
  # depth.bin.table
  # shallow    deep  step
  # [1,]    0.00  157.50  2.50
  # [2,]  157.50  161.25  3.75
  # [3,]  161.25  476.25  5.00
  # [4,]  476.25  483.75  7.50
  # [5,]  483.75 1113.75 10.00
  # [6,] 1113.75 1128.75 15.00
  # [7,] 1128.75 2388.75 20.00
  
  resid.bins <- c(0,60,120,180,10000)
  resid.table.colnames <- list("min", "max","R5")
  resid.table <- matrix(data=0, nrow=4, ncol=3, dimnames=list(NULL,resid.table.colnames))
  resid.table[,1] <- resid.bins[1:4]
  resid.table[,2] <- resid.bins[2:5]
  resid.table[,3] <- resid.bins[2:5]/10
  
  ################################################################################################
  # Prepare the abstracted or detailed time-depth data to be processed
  ################################################################################################
  
  # set number of columns in the ds matrix	
  ds.col <- breakpoints+1
  # create vector to hold faulty dives for which this code won't work. these are very short dives where breakpoints happen adjacent to each other and the slopes of segments cannot be worked out
  faulty <- NULL
  
  if (BS.dive==FALSE){
    
    if(missing(detailed.depth)){cat("Error: detailed.depth data are missing")}
    dive.time <- c(0:(length(detailed.depth)-1)) # time starts at 0
    dive.depth <- detailed.depth
    detailed.data <- as.data.frame(cbind(dive.time=dive.time, dive.depth=dive.depth))
    
    true.times <- detailed.data$dive.time*res/60 # res gives the resolution of the detailed dive data in seconds. the default is one time-depth pair every 4 seconds.
    true.depths <- detailed.data$dive.depth
    # depth stayed positive right until it gets plotted
    
  } else {
    
    if(missing(divetable)){cat("Error: divetable is missing")}	
    divedur.col <- grep("DIVE", names(divetable))
    
    tt <- as.numeric(divetable[dn,which(names(divetable)=="T1"):which(names(divetable)=="T4")]) # percent time where BS points occur
    bs.dive.time <- ((c(0, tt, 100)*0.01)*divetable[dn,divedur.col]) # these are the non-interpolated broken-stick time points
    bs.dive.timeIND <- round(bs.dive.time+1, digits=2) # these are the non-interpolated broken-stick time points turned into an index (i.e., starting from 1 and ending in dive.dur+1)
    
    dd <- as.numeric(divetable[dn,which(names(divetable)=="D1"):which(names(divetable)=="D4")])
    bs.dive.depth <- c(0, dd, 0) # these are the non-interpolated broken-stick depth points
    
    bs.round.times <- round(bs.dive.time, digits=2) # rounding is needed because this represents time since the dive started, for plotting. it is not an index
    bs.round.depths <- round(bs.dive.depth, digits=2)
    
    bin.row.index <- vector("numeric", length=breakpoints)
    step.size.vec <- vector("numeric", length=breakpoints)
    
    # work out the depth bin for each of the breakpoints. this is determined by the depth itself.
    for(aa in 2:(breakpoints+1)){ 
      bin.row.index[aa-1] <- which(bs.round.depths[aa] >= depth.bin.table[,1] & bs.round.depths[aa] < depth.bin.table[,2]) 
      step.size.vec[aa-1] <- depth.bin.table[bin.row.index[aa-1], 3]
    }
    step.size.vec <- c(2.5,step.size.vec,2.5)
    max.step.size <- max(step.size.vec)
    
    # work out the resolution at which the tag stored the dive
    if (divetable[dn,divedur.col] <= 256*4) {
      sampletime.steps <- floor(divetable[dn,divedur.col]/4)
      sampletime.res <- 4
    } else {
      if (divetable[dn,divedur.col] <= 256*8) {
        sampletime.steps <- floor(divetable[dn,divedur.col]/8)
        sampletime.res <- 8
      } else {
        if(divetable[dn,divedur.col] <= 256*16) {
          sampletime.steps <- floor(divetable[dn,divedur.col]/16)
          sampletime.res <- 16
        } else {
          sampletime.steps <- floor(divetable[dn,divedur.col]/32)
          sampletime.res <- 32
        }
      }
    }
    
    sample.t <- vector("numeric", length=1)
    seg.length <- vector("numeric", length=ds.col)
    # work out the time points at which to interpolate the dive depths at the appropriate resolution
    for (a in 1:(breakpoints+1)){ # a <- 5
      
      seg.length[a] <- round((bs.dive.timeIND[a+1]-bs.dive.timeIND[a])/round(sampletime.res))	
      if(seg.length[a]<1){seg.length[a] <- 1}
      
      seg.sample.times <- round(seq(from=bs.dive.timeIND[a], to=bs.dive.timeIND[a+1], by=(bs.dive.timeIND[a+1]-bs.dive.timeIND[a])/seg.length[a]), digits=3) 
      
      if (a==1){ # if you're working out times along the first segment you need to start it at 1 otherwise you end up with the wrong number of time points. you proceed to overwrite the value in the location of the last slot for each segment from here on
        sample.t[a:(a+seg.length[a])] <- seg.sample.times
      } else {
        
        sample.t[length(sample.t):(length(sample.t)+seg.length[a])] <- seg.sample.times
        
      }
      
    }
    
    sample.depths <- numeric(length(sample.t))
    for (e in 1:length(bs.dive.time)) { # e <-5
      sample.depths[round(sample.t,digits=2)==round(bs.dive.timeIND[e],digits=2)] <- bs.round.depths[e]
    }
    
    sample.times <- round(sample.t,digits=3) # work out the location of the depths on the original time scale (sec) 
    sample.index <- c(1:length(sample.t)) 
    
  } # this closes the else statement for if(BS.dive==FALSE)
  
  ################################################################################################
  # Apply the abstraction algorithm to my true or simulated dive: here using the Broken Stick Model
  ################################################################################################
  
  # create a matrix with named columns that holds the times and depths 
  if (BS.dive==FALSE){
    
    time.depth <- matrix(data=c(true.times, true.depths), nrow=length(true.times), ncol=2, dimnames=list(NULL,c("times", "depths"))) 
    times <- true.times
    depths <- true.depths
    
  } else {
    
    BS.time.depth <- matrix(data=c(sample.times, sample.depths), nrow=length(sample.times), ncol=2, dimnames=list(NULL,c("times", "depths"))) 
    times <- sample.times
    depths <- sample.depths
    
  } # this closes the else statement for if(BS.dive==FALSE) 
  
  # we need to iterate ii times so this will give us ii columns of residuals, where the last column, e.g., d5, is the residual at the next iteration that would have resulted in the addition of another point. ds.col is defined in the function arguments, depends on how many internal points you want the broken-stick to calculate
  ii <- ds.col 
  
  # this dataframe holds the residuals caculated at each iteration of the broken-stick algorithm. zeros represent maxima which correspond to selected depth points (i.e. D1 D2 D3 D4 etc) calculated for each segment
  ds.colnamestring <- paste("d", c(1:ds.col), sep="") 
  ds <- matrix(data=0, nrow=length(times), ncol=ii, dimnames=list(NULL,ds.colnamestring))
  
  R.name <- paste("R", c(1:ds.col), sep="") # creates a vector of residual names, R1:Rii
  # [1] "R1" "R2" "R3" "R4" "R5" etc
  R.vec <- matrix(0, nrow=ds.col, ncol=1, dimnames=list(R.name))
  
  # this 3D array holds the coordinates of depth points marking the start and end points for each segment as we break the "stick" in different places
  s <- array(0, dim=c(ii+1,ii+1,2))
  s[1,1,1] <- 1
  s[1,1,2] <- length(times) # this essentially gives your dive duration (DIVE.DUR or DIVE_DUR)
  
  if (BS.dive==TRUE){
    
    # If the dive data are abstracted and have been produced by the broken-stick algorithm then interpolate between broken-stick points to get a piecewise linear dive.
    # The dloc vector has to start with 1 not 0, because it is also an index. 
    dloc <- c(1, which(depths!=0), length(times)) # locations of bs.depths in vector of interpolated depths
    x <- times
    y <- vector(mode="numeric", length=length(times)) # length=max(times)+1)
    
    for (b in 1:(breakpoints+1)){ # because for 4 points you will get 5 segments
      
      up <- dloc[b+1]
      lo <- dloc[b]
      
      seg.slope <- round((depths[up]-depths[lo])/(times[up]-times[lo]), digits=3)
      # find the constant (b) in the equation of the line defined by the segment, if the b = y - mx, using y and x at the beginning of the segment 
      seg.b <- round((depths[lo] - seg.slope*times[lo]), digits=3)
      # create a numeric vector the length of index to hold the y values along the line segment joining the beginning and the end of the segment
      seg.ind <- c(lo:up)
      seg.x <- seq(from=times[lo], to=times[up], by=(times[up]-times[lo])/(up-lo)) # this will be 2 digits cos of the precision of sampletimes that I set at the top
      seg.y <- vector(mode="numeric", length=length(seg.x))
      
      # loop through the values in index and calculate a y for each one, along the line segment
      for (c in 1:length(seg.x)){ # 
        seg.y[which(seg.x==seg.x[c])] <- round((seg.slope*seg.x[c] + seg.b), digits=3)
        # abstracted depths along line segment
      } 				
      if(b==breakpoints+1){seg.y[length(seg.y)] <- 0}
      # round(sample.t,digits=2)==round(bs.dive.timeIND[e],digits=2) just for demo
      y[lo:up] <- seg.y
      
    } # closes for (b in...) loop
    
    y[y<0] <- 0
    BS.time.depth <- data.frame(times=times/60, depths=round(y, digits=2))	# THIS IS BEING RE-ASSIGNED HERE
    
  } else {
    
    x <- times
    y <- depths
    BS.time.depth <- data.frame(times=x, depths=y)	
    
  } # closes if (BS.dive==TRUE) else statement
  
  ############################################### BROKEN-STICK
  # this next part of the function is intended to emulate the broken-stick algorithm and produde an abstracted time-depth dive profile from a detailed profile that is assumed to be truth
  # the i loop runs through the number of iterations corresponding to the break-points and line segments you wish to produce (in our case 4 points, 5 segments - including dive start to dive end segment)
  # the j loop runs through the number of line segments produced at each iteration of i finding the distances from the line segment in question and true path, selects the maxima and then compares them. the biggest of those maximum departures is selected as the next break point. e.g. when i=2 j will go from 1 to 2 and find the distances from the true profile to segments drawn so far.
  #
  # T.Photopoulou, P. Lovell & J.Matthiopoulos, June 2014
  ###############################################
  
  for (i in 1:ii){ # where ii is equal to ds.col, so that i runs from 1 to 5, forming 5 segments
    # i <- 5
    dmax <- rep(0,i) # dmax is going to be a vector of maximum residuals (departures) for the segments being considered
    imax <- rep(0,i) # imax is going to be a vector of the row indeces for the maximum residuals (departures) for the segments being considered
    
    for (j in 1:i){
      # j runs from 1 to i. so in iteration 3 (i=3) j will run from 1 to 3 creating the indeces for the start/end of the line segments and calculating the max residual for each segment  
      index <- s[i,j,1]:s[i,j,2]
      # index only gives the row numbers of the (time) points across which segments are to be drawn (gives a sequence not just 2numbers), you also need the depths (y coords)
      temp <- V.distances(index, y[index]) # y is used here rather than depths, these are the interpolated depths along the line segments between breakpoints 
      dists <- temp$dists
      # dists is a vector of distances corresponding to the segment defined by index (x coords) and depths[index] (y coords)
      ds[index, i] <- dists
      # put the distances contained in dists into the appropriate column in the ds dataframe which is set up to store distances arising at each iteration of i 
      if(!length(index[which.max(dists)])){faulty <- c(faulty,i); next(paste("Sorry, dive", dn, "is faulty"))}
      imax[j] <- index[which.max(dists)]
      # imax is a vector containing the x coordinate for the maximum distance
      dmax[j] <- max(dists)
      # dmax is a vector containing the actual value for the maximum distance (for i=1,j=1 this is in fact the depth in metres)
      
      # now I have two loops (one for i (ii iterations), one for j within each i (runs from 1 to i, so that it checks all segments each time)) and one function that works out distances for all segments defined within each j
      # next the biggest dmax needs to be picked out and matched to the corresponding time point, then the break needs to be made, and the segment coordinates added to the 3D array called s
      
    } # this bracket closes the for (j in 1:i) loop, i is still in play 
    
    k <- which.max(dmax)
    # k is the position (element number) of the maximum distances in the dmax vector which has length i
    # when i=1 k is always going to be 1 as there will only be one element to dmax, when i=2 k could be 1 or 2, when i=3 k could be 1, 2 or 3 etc.
    ibreak <- imax[k]
    # ibreak is the biggest residual value among the maxima held in the vector imax which has length i 
    
    v <- s[i,,]
    # this is a dummy variable that contains the first row of the s array, e.g. a 2D array of dim(5,2) where 5 is the y dim and 2 is the z dim of the 3D array s
    vstart <- v[k,1]
    # vstart is the beginning of segment where the max occurs, since v contains the index for the segment's start and end points
    vend <- v[k,2]
    # vend is the end of segment where the max occurs, since v contains the index for the segment's start and end points
    v[k,] <- c(vstart, ibreak)
    # replace the two elements in k with the original vstart and the breakpoint to form the first of two new segments
    v[i+1,] <- c(ibreak, vend)
    # make a new entry in the next open position in the ith row of s (i+1) starting with the breakpoint (ibreak) and ending with the original vend to form the second of two new segments
    # this is why you need ii+1 spaces in the s array. 
    s[i+1,,] <- v
    # place v in the next row down (of the s array)
    R.vec[i] <- max(ds[,i])
    
  } # this bracket closes the for(i in 1:ii) loop
  
  BS.order <- c(0,order(s[ii,2:ii,1]),ds.col) # this is the order in which the BSM points were added to the profile if you are looking at a plot of the time-depth profile. it will read something like 0 4 1 3 2 5. this means that outside the first and last points (always 0 and 5, respectively) the second breakpoint was added to the profile first. 
  
  if (BS.dive==TRUE){
    time.depth <- NULL
    if(breakpoints==4 && ds.col==5){
      R5.row.index <- which(divetable$RESIDUAL[dn] >= resid.table[,1] & divetable$RESIDUAL[dn] < resid.table[,2])
      R5.trans <- resid.table[R5.row.index, 3]
      
      R.vec[ii] <- R5.trans 
      
    } # closes the if(breakpoints...) statement
    
  } else { # this closes the if(BS.dive=TRUE) statement
    
    sample.index <- seq(1:nrow(time.depth))
    R5.trans <- dn <- divetable <- step.size.vec <- NULL
    sampletime.res <- res
    dloc <- c(1, s[ds.col,c(1:ds.col),2]) # should this perhaps be c(s[ds.col,c(1:ds.col),1], max(dive.time))
    
    sample.t <- vector("numeric", length=1)
    seg.length <- vector("numeric", length=ds.col)
    sample.depths <- numeric(length(sample.t))
    s.dloc <- sort(dloc)
    
    # this section generates interpolated depths along the abstracted dive segments
    for (a in 1:(breakpoints+1)){ # a <- 5
      
      seg.length[a] <- length(s.dloc[a]:s.dloc[a+1])
      if(seg.length[a]<1){seg.length[a] <- 1}
      
      seg.sample.times <- seq(from=s.dloc[a], to=s.dloc[a+1], by=1) 
      seg.sample.depths <- seq(from=depths[s.dloc[a]], to=depths[s.dloc[a+1]], length.out=seg.length[a]) 
      
      if (a==1){ # if you're working out times along the first segment you need to start it at 1 otherwise you end up with the wrong number of time points. you proceed to overwrite the value in the location of the last slot for each segment from here on. 
        sample.t[a:seg.length[a]] <- seg.sample.times
      } else {
        
        sample.t[seg.sample.times[1]:(max(seg.sample.times))] <- seg.sample.times
        
      }
      
      sample.depths[seg.sample.times[1]:(max(seg.sample.times))] <- seg.sample.depths		
      
    } # closes a loop
    
    BS.time.depth <- data.frame(times=true.times, depths=sample.depths) # this is being assigned for the first time here for detailed dives  
    
    bin.row.index <- vector("numeric", length=breakpoints)
    step.size.vec <- vector("numeric", length=breakpoints)
    for(aa in 2:(breakpoints+1)){ 
      bin.row.index[aa-1] <- which(BS.time.depth[s[aa,aa,1],2] >= depth.bin.table[,1] & BS.time.depth[s[aa,aa,1],2] < depth.bin.table[,2]) # still not sure if this should be && or not....
      step.size.vec[aa-1] <- depth.bin.table[bin.row.index[aa-1], 3]
    } # closes aa loop
    step.size.vec <- c(2.5,step.size.vec,2.5)
    max.step.size <- max(step.size.vec)
    
  } # this closes the if(BS.dive=TRUE) "else" statement
  
  ############################################### PLOTTING
  
  if (plot==TRUE){
    
    if (dev.new==TRUE){
      dev.new()
    }
    
    if (plot.truth==TRUE){	
      
      plot(true.times, -true.depths, ylim=c(-(max(true.depths)+0.1*(max(true.depths))), 0.1*(max(true.depths))), ylab="Depth (m)", xlab="Time (min)", col=2, type="l", cex.axis=cex.axis, lty=2, cex.lab=cex.lab)
    } # closes plot.truth=TRUE 
    
    if (plot.bsm==TRUE){
      if(plot.truth==FALSE){
        plot(BS.time.depth$times, -BS.time.depth$depths, ylim=c(-(max(BS.time.depth$depths)+0.1*(max(BS.time.depth$depths))), 0.1*(max(BS.time.depth$depths))), ylab="Depth (m)", xlab="Time (min)", col=2, type="n", cex.axis=cex.axis, lty=2, cex.lab=cex.lab) 
        
      } 
      
      for (i in 1:ds.col){
        segments(x0=BS.time.depth[s[ds.col,i,1],1], y0=-BS.time.depth[s[ds.col,i,1],2], x1=BS.time.depth[s[ds.col,i,2],1], y1=-BS.time.depth[s[ds.col,i,2],2])
      }
      
      if (draw.numbers == TRUE){
        text(x=BS.time.depth[s[1,1,1],1],y=-(BS.time.depth[s[1,1,1],2]-0.1*(max(BS.time.depth$depths))),"0", cex=cex.lab) 
        text(x=BS.time.depth[s[1,1,2],1],y=-(BS.time.depth[s[1,1,2],2]-0.1*(max(BS.time.depth$depths))),paste(ds.col), cex=cex.lab)
        
        for (i in 2:(ds.col)){
          text(x= BS.time.depth[s[i,i,1],1],y=-(BS.time.depth[s[i,i,1],2]+0.1*(max(BS.time.depth$depths))),paste(i-1), cex=cex.lab) 
        }
        
      } # this bracket closes the if(draw.numbers==TRUE) statement
      
      if (draw.lines == TRUE){
        abline(v=BS.time.depth[which.max(BS.time.depth$depths),1], h=c(0, -max(BS.time.depth$depths)), col=8)
      } # this bracket closes the if(draw.lines==TRUE) statement
      
    } # this bracket closes the if(plot.bsm==TRUE) statements
    
  } # this bracket closes the if(plot==TRUE) statement
  
  out <- list(ds=ds, s=s, time.depth=time.depth, BS.time.depth=BS.time.depth, seg.length=seg.length, R.vec=R.vec, R5.trans=R5.trans, dive.number=dn, sample.index=sample.index, breakpoints=breakpoints, sampletime.res=sampletime.res, step.size.vec=step.size.vec, BS.order=BS.order, dloc=dloc, divetable=divetable, dn=dn, depth.bin.table=depth.bin.table)
  return(out)
  
}

################ find.Xzone
# This function calculates the dive zone for a dive that has been abstracted using the broken-stick algorithm 

find.Xzone <- function(output, 
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
      
      for(i in 1:(R.index+1)){
        
        segments(x0=BS.time.depth[s[R.index+1,i,1],1], y0=-BS.time.depth[s[R.index+1,i,1],2], x1=BS.time.depth[s[R.index+1,i,2],1], y1=-BS.time.depth[s[R.index+1,i,2],2], lwd=lwd-1)
        
      }
      
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
  
  out <- list(tt=tt, allX.lines=allX.lines, BS.time.depth=BS.time.depth, Xzones=Xzones, R.index=R.index, Resid.index=Resid.index, R.vec=R.vec, dive.number=dn, DZI=dz.index) 
  return(out)
  
}

################ xzone.wrapper
# This wrapper function lets you run the abstraction and/or the dive zone function(s) on whole datasets of dives, in "SMRU divetable" form

xzone.wrapper <- function(BS.dive=TRUE, divetable=divetable, start.dive, end.dive, breakpoints=4, plot.bsm=FALSE, plot.xzone=TRUE, cex.axis=1.3, cex.lab=1.3){
  
  ## ARGUMENTS TO FUNCTION - use args(xzone.wrapper) to see the defaults
  
  # BS.dive: the object holding the output from get.BSMdive 
  # divetable: is a logical argument indicating whether the dive data have already been abstracted by the broken-stick model
  # start.dive: is a logical argument and indicates whether you want your dive to be plotted
  # end.dive: is a logical argument and indicates whether you want ta legend to be plotted
  # breakpoints: is an integer and is the number of breakpoints you want in your abstracted dive. the default is 4
  # plot.bsm: is a logical argument and indicates whether you want the broken-stick abstracted dive to be plotted
  # plot.xzone: is a logical argument and indicates whether you want the dive zone to be plotted
  # cex.axis: is a real number and indicates the size of the numbers in the axes
  # cex.lab: is a real number and indicates the size of the axis labels 
  
  diveseq <- start.dive:end.dive	
  faulty <- vector("logical", length(diveseq))	
  Rt.vec <- dzi <- vector("numeric", length(diveseq))	 
  divetimes <- divedepths <- diveid <- up <- lo <- vector("list", length=length(diveseq))
  starttime <- maxdep.vec <- tsamples.vec <- vector("numeric", length=length(diveseq))
  xlines <- vector("list", length=length(diveseq))
  
  for (i in 1:length(diveseq)){
    
    cat("Dive", i, "of", length(diveseq), "\n")
    flush.console()
    
    bsm <- try(get.BSMdive(BS.dive=TRUE, divetable=divetable, breakpoints=breakpoints, dn=diveseq[i], dev.new=FALSE, plot=plot.bsm, plot.truth=FALSE, cex.axis=cex.axis, cex.lab=cex.lab), silent=TRUE)
    
    if (class(bsm)=="try-error") {faulty[i] <- TRUE; dzi[i] <- "NA"; Rt.vec[i] <- "NA"} else {
      
      divetimes[[i]] <- bsm$BS.time.depth$times
      divedepths[[i]] <- bsm$BS.time.depth$depths
      diveid[[i]] <- rep(paste("dive",diveseq[i],sep=""), nrow(bsm$BS.time.depth))
      starttime[i] <- as.POSIXct(divetable$DE.DATE[i]-divetable$DIVE.DUR[i])
      ts <- as.POSIXct(starttime[i], origin="1970-01-01 00:00:00")+divetimes[[i]]*60
      divetimes[[i]] <- ts
      maxdep.vec[i] <- bsm$divetable$MAX.DEP[i]
      tsamples.vec[i] <- length(bsm$BS.time.depth$times)
      
      xzone <- find.Xzone(output=bsm, BS.dive=BS.dive, plot=plot.xzone, dev.new=FALSE, plot.truth=FALSE)
      up[[i]] <- -xzone$Xzones[,(xzone$R.index*2)-1]
      lo[[i]] <- -xzone$Xzones[, xzone$R.index*2]
      xlines[[i]] <- matrix(data=cbind(up[[i]], lo[[i]]), ncol=2, dimnames=list(NULL, c("dz.UP", "dz.LO")))
      
      dzi[i] <- as.numeric(xzone$DZI)
      Rt.vec[i] <- as.numeric(bsm$R.vec[breakpoints,])
      
    }
    
  }	
  names(xlines) <- paste("dive",diveseq,sep="")
  out.t <- data.frame(ref=divetable$ref[diveseq], dive=diveseq, max.dep=maxdep.vec, n.samples=tsamples.vec, Rmax=as.numeric(Rt.vec), DZI=as.numeric(dzi))
  out <- list(dives=out.t, divezone.limits=xlines, faulty=which(faulty))
  return(out)
  
}

############################################### VERTICAL DISTANCES
# This function calculates vertical distances between two time-series of depth points, 
# one linear and one non-linear. 
# It takes two arguments, index and dep. The first argument, index, is a vector including the sequence of 
# time (x) values, and the second argument, dep, is a vector of the depth (y) values that 
# correspond to the time points in index. Using the slope and end points of the linear time-series
# to calculate points along that line (y = mx + b) so that I can take the difference between 
# the line and the true dive path at each time point i. 
#
# T.Photopoulou, March 2010
###############################################

V.distances <- function(index, dep){
  
  # find the slope of the segments, m=(y2-y1)/(x2-x1)
  seg.slope <- (dep[length(index)]-dep[1])/(max(index)-min(index))
  # find the constant (b) in the equation of the line defined by the segment, if the b = y - mx, using y and x at the beginning of the segment 
  seg.b <- dep[1] - seg.slope*min(index)
  # create a numeric vector the length of index to hold the y values along the line segment joining the beginning and the end of the segment
  seg.y <- vector(mode="numeric", length=length(index))
  
  # loop through the values in index and calculate a y for each one, along the line segment
  for (i in index){
    seg.y[which(index==i)] <- seg.slope*i + seg.b
    # abstracted depths along line segment
  }	
  
  # find the absolute difference between the observed depth and the line segment
  temp.dists <- abs(dep-seg.y)
  out <- list(dists=temp.dists, seg.y=seg.y)
  return(out)  		
}

