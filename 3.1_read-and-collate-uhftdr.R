
#### Include dive data per time step as a propotion of water column ####

## 2016 ##

filenames <- list.files(path="Data/Depth Data/Files/2016", pattern="*.tdr", full.names=TRUE)

for(i in 1:length(filenames)) {
  X <- read.csv(filenames[i], header=FALSE, sep="") # header=FALSE creates 14 separate columns for each variable rather than one column for everything 
  skip = 3 #skips first 3 lines of data so you're just left with two rows as headers which are full of NAs and therefore easily subsetable 
  
  X <- na.omit(X) # removes the first 2 rows with NA values so you're just left with the data without headers
  
  colnames(X) <- c("Dive Number","Year", "Month","Day", "Hour", "Min","Sec", "Depth","Dummy", "Temperature") # adds names of columns
  
  ds<-paste(filenames[i],sep="") #adds the filename to ds
  
  ds<-substr(ds, 45, nchar(ds)-4) # remove the last 13 charaters (-COMPLETE.pos)
  
  dsn <- paste(ds, i, sep=".")
  
  X$tag<- substr(ds, 4, 8) # removes "Tag" from the tag column entry (begin at the cut at the 4th character and end at the 8th character)
  
  assign(dsn, X) # Give the name "ds" a value of X (data.frame)
  
  print(i)
  
}


filenames <- list.files(path="Data/Depth Data/Files/2017", pattern="*.tdr", full.names=TRUE)



for(i in 1:length(filenames)) {
  X <- read.csv(filenames[i], header=FALSE, sep="") # header=FALSE creates 14 separate columns for each variable rather than one column for everything 
  # skip = 3 skips first 3 lines of data so you're just left with two rows as headers which are full of NAs and therefore easily subsetable 
  
  X <- na.omit(X) # removes the first 2 rows with NA values so you're just left with the data without headers
  
  colnames(X) <- c("Year", "Month","Day", "Hour", "Min","Sec", "Depth","Dummy", "Temperature", "Voltage") # adds names of columns
  
  ds<-paste(filenames[i],sep="") #adds the filename to ds
  
  ds<-substr(ds, 45, nchar(ds)-4) # remove the last 13 charaters (-COMPLETE.pos)
  
  dsn <- paste(ds, i, sep=".")
  
  X$tag<- substr(ds, 4, 8)
  
  assign(dsn, X) # Give the name "ds" a value of X (data.frame)
  
  print(i)
  
}


filenames <- list.files(path="Data/Depth Data/Files/2018", pattern="*.tdr", full.names=TRUE)



for(i in 1:length(filenames)) {
  X <- read.csv(filenames[i], header=FALSE, sep="") # header=FALSE creates 14 separate columns for each variable rather than one column for everything 
  # skip = 3 skips first 3 lines of data so you're just left with two rows as headers which are full of NAs and therefore easily subsetable 
  
  X <- na.omit(X) # removes the first 2 rows with NA values so you're just left with the data without headers
  
  colnames(X) <- c("Year", "Month","Day", "Hour", "Min","Sec", "Depth","Dummy", "Temperature", "Voltage") # adds names of columns
  
  ds<-paste(filenames[i],sep="") #adds the filename to ds
  
  ds<-substr(ds, 45, nchar(ds)-4) # remove the last 13 charaters (-COMPLETE.pos)
  
  dsn <- paste(ds, i, sep=".")
  
  X$tag<- substr(ds, 4, 8) # removes "Tag" from the tag column entry (begin at the cut at the 4th character and end at the 8th character)
  
  assign(dsn, X) # Give the name "ds" a value of X (data.frame)
  
  print(filenames[i])
  
}


## Multiple files for a single tag so need to combine to form single data frame for each tag:

tdr51009 <- rbindlist(mget(ls(pattern = "51009"))) 
tdr51011 <- rbindlist(mget(ls(pattern = "51011")))
tdr51019 <- rbindlist(mget(ls(pattern = "51019")))
tdr51020 <- rbindlist(mget(ls(pattern = "51020")))
tdr51022 <- rbindlist(mget(ls(pattern = "51022")))
tdr51025 <- rbindlist(mget(ls(pattern = "51025")))
tdr51026 <- rbindlist(mget(ls(pattern = "51026")))
tdr51029 <- rbindlist(mget(ls(pattern = "51029")))
tdr51030 <- rbindlist(mget(ls(pattern = "51030")))
tdr51031 <- rbindlist(mget(ls(pattern = "51031")))

tdr51104 <- rbindlist(mget(ls(pattern = "51104")))
tdr51120 <- rbindlist(mget(ls(pattern = "51120")))
tdr51105 <- rbindlist(mget(ls(pattern = "51105")))
tdr51109 <- rbindlist(mget(ls(pattern = "51109")))
tdr51111 <- rbindlist(mget(ls(pattern = "51111")))
tdr51100 <- rbindlist(mget(ls(pattern = "51100")))
tdr51101 <- rbindlist(mget(ls(pattern = "51101")))
tdr51112 <- rbindlist(mget(ls(pattern = "51112")))
tdr51119 <- rbindlist(mget(ls(pattern = "51119")))
tdr51114 <- rbindlist(mget(ls(pattern = "51114")))
tdr51116 <- rbindlist(mget(ls(pattern = "51116")))
tdr51115 <- rbindlist(mget(ls(pattern = "51115")))
tdr51108 <- rbindlist(mget(ls(pattern = "51108")))
tdr51117 <- rbindlist(mget(ls(pattern = "51117")))


tdr51102 <- rbindlist(mget(ls(pattern = "51102")))
tdr51110 <- rbindlist(mget(ls(pattern = "51110")))
tdr51118 <- rbindlist(mget(ls(pattern = "51118")))
tdr51121 <- rbindlist(mget(ls(pattern = "51121")))
tdr51122 <- rbindlist(mget(ls(pattern = "51122")))
tdr51124 <- rbindlist(mget(ls(pattern = "51124")))
tdr51125 <- rbindlist(mget(ls(pattern = "51125")))
tdr51126 <- rbindlist(mget(ls(pattern = "51126")))
tdr51127 <- rbindlist(mget(ls(pattern = "51127")))
tdr51128 <- rbindlist(mget(ls(pattern = "51128")))
tdr51129 <- rbindlist(mget(ls(pattern = "51129")))
tdr51130 <- rbindlist(mget(ls(pattern = "51130")))
tdr51131 <- rbindlist(mget(ls(pattern = "51131")))
tdr51132 <- rbindlist(mget(ls(pattern = "51132")))
tdr51134 <- rbindlist(mget(ls(pattern = "51134")))
tdr51136 <- rbindlist(mget(ls(pattern = "51136")))



TDR2016 <- rbindlist(mget(ls(pattern = "tdr510")))
TDR2017 <- rbind(tdr51104,tdr51120,tdr51105,tdr51109,tdr51111,tdr51100,tdr51101,tdr51112,tdr51119,
                 tdr51114,tdr51116,tdr51115,tdr51108,tdr51117) 
TDR2018 <- rbind(tdr51102,tdr51110,tdr51118,tdr51121,tdr51122,tdr51124,tdr51125,
                 tdr51126,tdr51127,tdr51128,tdr51129,tdr51130,tdr51131,tdr51132,tdr51134,tdr51136)



TDR2016$datetime<- paste(TDR2016$Year,"/",TDR2016$Month,"/",TDR2016$Day," ",TDR2016$Hour,":",TDR2016$Min,":",TDR2016$Sec, sep="")
temp <- as.POSIXct(as.character(TDR2016$datetime), "%Y/%m/%d %H:%M:%S", tz="UTC")
TDR2016$datetime <- temp
temp2 <- as.numeric(difftime(temp, "2016/09/28 12:00:00", tz="UTC", units="secs"))
TDR2016$secs<- temp2
rm(temp, temp2)


TDR2016 <- subset(TDR2016, Depth>-0.01 & Depth<200 & Dummy < 200)
TDR2016$Month <- as.Date(TDR2016$datetime)


TDR2017$datetime<- paste(TDR2017$Year,"/",TDR2017$Month,"/",TDR2017$Day," ",TDR2017$Hour,":",TDR2017$Min,":",TDR2017$Sec, sep="")
temp <- as.POSIXct(as.character(TDR2017$datetime), "%Y/%m/%d %H:%M:%S", tz="UTC")
TDR2017$datetime <- temp
temp2 <- as.numeric(difftime(temp, "2016/09/28 12:00:00", tz="UTC", units="secs"))
TDR2017$secs<- temp2
rm(temp, temp2)


# Clean data to include only positive depth and logical depth values: Dummy <200 is same reliance metric as residual error
TDR2017 <- subset(TDR2017, Depth>-0.01 & Depth<200 & Dummy < 200)
TDR2017$Month <- as.Date(TDR2017$datetime)

TDR2018$datetime<- paste(TDR2018$Year,"/",TDR2018$Month,"/",TDR2018$Day," ",TDR2018$Hour,":",TDR2018$Min,":",
                         TDR2018$Sec, sep="")
temp <- as.POSIXct(as.character(TDR2018$datetime), "%Y/%m/%d %H:%M:%S", tz="UTC")
TDR2018$datetime <- temp
temp2 <- as.numeric(difftime(temp, "2016/09/28 12:00:00", tz="UTC", units="secs"))
TDR2018$secs<- temp2
rm(temp, temp2)


# Clean data to include only positive depth and logical depth values: Dummy <200 is same reliance metric as residual error
TDR2018 <- subset(TDR2018, Depth>-0.01 & Depth<200 & Dummy < 200)
TDR2018$Month <- as.Date(TDR2018$datetime)


save(TDR2016, file="Data/Depth Data/TDR2016.Rd")
save(TDR2017, file=paste("Data/Depth Data/TDR2017.Rd"))
save(TDR2018, file=paste("Data/Depth Data/TDR2018.Rd"))

