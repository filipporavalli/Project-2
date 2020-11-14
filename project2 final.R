
#Please run all code up to and including "Final Algorithm" for testing predictions
  ##In regression tab, only necessary to run "week1long", "week1lat", "week1elapsedtime"

##### Load libraries ####
library(data.table)
library(jsonlite)
library(sp)
library(rgdal)
library(dlm)
library(plyr)
library(ggplot2)
library(ggpubr)
library(dplyr)
### end


##### Read data in ####
data18<-read_json("~/Desktop/gps/20200818114606.geojson")
data19<-read_json("~/Desktop/gps/20200819132607.geojson")
data20<-read_json("~/Desktop/gps/20200820151044.geojson")
data21<-read_json("~/Desktop/gps/20200821111447.geojson")
data24<-read_json("~/Desktop/gps/20200824130857.geojson")
data25<-read_json("~/Desktop/gps/20200825121346.geojson")
data26<-read_json("~/Desktop/gps/20200826131614.geojson")
data27<-read_json("~/Desktop/gps/20200827113234.geojson")
data28<-read_json("~/Desktop/gps/20200828122627.geojson")
data29<-read_json("~/Desktop/gps/20200828130816.geojson")
data30<-read_json("~/Desktop/gps/20200831115147.geojson")
### end


##### Combine Data ####
geodf<-data.frame(long=NA, lat=NA, time_stamp=NA)
geodf1<-data.frame()
i<-1
for (i in 1:length(data18[["features"]])){
  geodf$long[1]<-data18[["features"]][[i]][["geometry"]][["coordinates"]][[1]]
  geodf$lat[1]<-data18[["features"]][[i]][["geometry"]][["coordinates"]][[2]]
  geodf$time_stamp[1]<-data18[["features"]][[i]][["properties"]][["time"]]
  geodf1<-rbind(geodf1, geodf)
}

geodf<-data.frame(long=NA, lat=NA, time_stamp=NA)
geodf2<-data.frame()
i<-1
for (i in 1:length(data19[["features"]])){
  geodf$long[1]<-data19[["features"]][[i]][["geometry"]][["coordinates"]][[1]]
  geodf$lat[1]<-data19[["features"]][[i]][["geometry"]][["coordinates"]][[2]]
  geodf$time_stamp[1]<-data19[["features"]][[i]][["properties"]][["time"]]
  geodf2<-rbind(geodf2, geodf)
}

geodf<-data.frame(long=NA, lat=NA, time_stamp=NA)
geodf3<-data.frame()
i<-1
for (i in 1:length(data20[["features"]])){
  geodf$long[1]<-data20[["features"]][[i]][["geometry"]][["coordinates"]][[1]]
  geodf$lat[1]<-data20[["features"]][[i]][["geometry"]][["coordinates"]][[2]]
  geodf$time_stamp[1]<-data20[["features"]][[i]][["properties"]][["time"]]
  geodf3<-rbind(geodf3, geodf)
}


geodf<-data.frame(long=NA, lat=NA, time_stamp=NA)
geodf4<-data.frame()
i<-1
for (i in 1:length(data21[["features"]])){
  geodf$long[1]<-data21[["features"]][[i]][["geometry"]][["coordinates"]][[1]]
  geodf$lat[1]<-data21[["features"]][[i]][["geometry"]][["coordinates"]][[2]]
  geodf$time_stamp[1]<-data21[["features"]][[i]][["properties"]][["time"]]
  geodf4<-rbind(geodf4, geodf)
}

geodf<-data.frame(long=NA, lat=NA, time_stamp=NA)
geodf5<-data.frame()
i<-1
for (i in 1:length(data24[["features"]])){
  geodf$long[1]<-data24[["features"]][[i]][["geometry"]][["coordinates"]][[1]]
  geodf$lat[1]<-data24[["features"]][[i]][["geometry"]][["coordinates"]][[2]]
  geodf$time_stamp[1]<-data24[["features"]][[i]][["properties"]][["time"]]
  geodf5<-rbind(geodf5, geodf)
}

geodf<-data.frame(long=NA, lat=NA, time_stamp=NA)
geodf6<-data.frame()
i<-1
for (i in 1:length(data25[["features"]])){
  geodf$long[1]<-data25[["features"]][[i]][["geometry"]][["coordinates"]][[1]]
  geodf$lat[1]<-data25[["features"]][[i]][["geometry"]][["coordinates"]][[2]]
  geodf$time_stamp[1]<-data25[["features"]][[i]][["properties"]][["time"]]
  geodf6<-rbind(geodf6, geodf)
}

geodf<-data.frame(long=NA, lat=NA, time_stamp=NA)
geodf7<-data.frame()
i<-1
for (i in 1:length(data26[["features"]])){
  geodf$long[1]<-data26[["features"]][[i]][["geometry"]][["coordinates"]][[1]]
  geodf$lat[1]<-data26[["features"]][[i]][["geometry"]][["coordinates"]][[2]]
  geodf$time_stamp[1]<-data26[["features"]][[i]][["properties"]][["time"]]
  geodf7<-rbind(geodf7, geodf)
}

geodf<-data.frame(long=NA, lat=NA, time_stamp=NA)
geodf8<-data.frame()
i<-1
for (i in 1:length(data27[["features"]])){
  geodf$long[1]<-data27[["features"]][[i]][["geometry"]][["coordinates"]][[1]]
  geodf$lat[1]<-data27[["features"]][[i]][["geometry"]][["coordinates"]][[2]]
  geodf$time_stamp[1]<-data27[["features"]][[i]][["properties"]][["time"]]
  geodf8<-rbind(geodf8, geodf)
}


geodf<-data.frame(long=NA, lat=NA, time_stamp=NA)
geodf9<-data.frame()
i<-1
for (i in 1:length(data28[["features"]])){
  geodf$long[1]<-data28[["features"]][[i]][["geometry"]][["coordinates"]][[1]]
  geodf$lat[1]<-data28[["features"]][[i]][["geometry"]][["coordinates"]][[2]]
  geodf$time_stamp[1]<-data28[["features"]][[i]][["properties"]][["time"]]
  geodf9<-rbind(geodf9, geodf)
}

geodf<-data.frame(long=NA, lat=NA, time_stamp=NA)
geodf10<-data.frame()
i<-1
for (i in 1:length(data29[["features"]])){
  geodf$long[1]<-data29[["features"]][[i]][["geometry"]][["coordinates"]][[1]]
  geodf$lat[1]<-data29[["features"]][[i]][["geometry"]][["coordinates"]][[2]]
  geodf$time_stamp[1]<-data29[["features"]][[i]][["properties"]][["time"]]
  geodf10<-rbind(geodf10, geodf)
}


geodf<-data.frame(long=NA, lat=NA, time_stamp=NA)
geodf11<-data.frame()
i<-1
for (i in 1:length(data30[["features"]])){
  geodf$long[1]<-data30[["features"]][[i]][["geometry"]][["coordinates"]][[1]]
  geodf$lat[1]<-data30[["features"]][[i]][["geometry"]][["coordinates"]][[2]]
  geodf$time_stamp[1]<-data30[["features"]][[i]][["properties"]][["time"]]
  geodf11<-rbind(geodf11, geodf)
}

datafull<-rbind(geodf1, geodf2, geodf3, geodf4, geodf5, geodf6, geodf7, geodf8, geodf9, geodf10, geodf11)

### end

##### Fix time/location ####
datafull$time_stamp<-as.POSIXct(datafull$time_stamp, format="%Y-%m-%dT%H:%M:%S")


spat_df <- SpatialPointsDataFrame(coords=datafull[, c("long", "lat")],
                                  data=datafull['time_stamp'], 
                                  proj4string=CRS("+proj=longlat +datum=WGS84"))
# This step converts the longitude/latitude -> UTM
utm_df <- spTransform(spat_df, CRSobj = "+proj=utm +zone=12 +datum=WGS84")
utm_coords <- coordinates(utm_df)

#create dataframe w utm 
utmdatafull<-c()
utmdatafull$long<-utm_coords[,1]
utmdatafull$lat<-utm_coords[,2]
utmdatafull<-as.data.frame(utmdatafull)
utmdatafull$time_stamp<-utm_df@data[["time_stamp"]]



#Create elapsed time
utmdatafull<-utmdatafull[order(utmdatafull$time_stamp),]
i<-1
for (i in 1:nrow(utmdatafull)){
  if ((i-1)==0){
    utmdatafull$elapsedtime[i]<-NA
  }
  else{
    utmdatafull$elapsedtime[i]<- difftime(utmdatafull$time_stamp[i], (utmdatafull$time_stamp[i-1]), units="sec")
  }
}

utmdatafull$session[1]<-1
x<-1
i<-1
for (i in 2:nrow(datafull)) {
  if (utmdatafull$elapsedtime[i] > 120){
    x<-x+1
    utmdatafull$session[i]<-x
  }
  else {
    utmdatafull$session[i]<-x
  }}

##create net elapsed time
utmdatafull$elapsedtime[1]<-0
utmdatafull$totelapsed[1]<-0

i<-1
for (i in 2:nrow(utmdatafull)) {
  if (utmdatafull$elapsedtime[i]>120) {
    utmdatafull$totelapsed[i]<-0
  }
  else{
    utmdatafull$totelapsed[i]<-(utmdatafull$elapsedtime[i])+(utmdatafull$totelapsed[i-1])
  }
}

#Create subdata frames
dfweek1<-utmdatafull[which(mday(utmdatafull$time_stamp)<=22),]
dfweek2<-utmdatafull[which(mday(utmdatafull$time_stamp)>22),]

#Remove outlier sessions
utmtrain<-utmdatafull[!utmdatafull$session==22&!utmdatafull$session==15&
                        !utmdatafull$session==2&!utmdatafull$session==23,]

#Setting initial dataframe
intial<-utmtrain[utmtrain$totelapsed==0,]
i<-1
x<-1
for (i in 1:20){
  for(x in 1:unique(utmtrain$session)[i]){
    intial$medtime[i]<-median(utmtrain$totelapsed[utmtrain$session==x]) 
  }}
### end

##### Regression ####
#to/from
todata<-utmtrain[which(utmtrain$session==1|utmtrain$session==3|utmtrain$session==5|
                         utmtrain$session==7|utmtrain$session==9
                       |utmtrain$session==11 |utmtrain$session==13
                       |utmtrain$session==16 |utmtrain$session==18
                       |utmtrain$session==21),]

fromdata<-utmtrain[which(utmtrain$session==4|utmtrain$session==6|utmtrain$session==8|
                           utmtrain$session==10 |utmtrain$session==12 |utmtrain$session==14
                         |utmtrain$session==17 |utmtrain$session==20
                         |utmtrain$session==22 |utmtrain$session==24),]

ggplot(NULL, aes(long, lat))+
  geom_point(data=todata, color="red")+
  geom_point(data=fromdata, color="green")

mean(intial$medtime)
sd(finalduration)
#Regressions
week1long<-lm(log(long)~log(totelapsed)+totelapsed, 
              data=utmtrain[!utmtrain$totelapsed<913&!utmtrain$totelapsed>1697,])

week1lat<-lm(log(lat)~log(totelapsed)+long+totelapsed, 
             data=utmtrain[!utmtrain$totelapsed<913&!utmtrain$totelapsed>1697,])

week1elapsedtime<-lm(log(medtime)~log(lat)+log(long)+lat+long, 
                     data=intial[!intial$medtime<913 &!intial$medtime>1697,])

par(mfrow = c(2, 2))
plot(week1elapsedtime)

#testing
predictionlocation<-utmtrain[utmtrain$session==10&utmtrain$totelapsed==0,]

timepredict<-predict.lm(week1elapsedtime, predictionlocation, se.fit = TRUE, interval="prediction", level=0.90)
totelapsed<-exp(timepredict$fit[1])
timeprediction<-as.data.frame(totelapsed)

longpredict<-predict.lm(week1long, timeprediction, se.fit = TRUE, interval = "prediction", level = 0.90)
finallong<-data.frame(long=exp(longpredict$fit[1]), totelapsed=totelapsed)
latpredict<-predict.lm(week1lat, finallong, se.fit = TRUE, interval = "prediction", level= 0.90)
finallat<-as.data.frame(exp(latpredict$fit))
finalloc<-data.frame(long=finallong[,1], lat=finallat[,1], totelapsed=totelapsed)

predictintlong<-data.frame(upr=exp(longpredict$fit)[3], lwr=exp(longpredict$fit)[2])
predictintlat<-data.frame(upr=exp(latpredict$fit)[3], lwr=exp(latpredict$fit)[2])
predictinttime<-data.frame(upr=exp(timepredict$fit)[3], lwr=exp(timepredict$fit)[2])
class(predictintlat)
class(predictinttime)
subset(utmtrain, long>(finalloc$long-5) & long<(finalloc$long+5) & lat>(finalloc$lat-5) & lat<(finalloc$lat+5))

subset(utmtrain[utmtrain$session==10,], totelapsed>(finalloc$totelapsed-10) & totelapsed<(finalloc$totelapsed+10))
### end

##### Final Algorithm ####
testgps<-read_json("~/Desktop/test_gps/20200901112100.geojson")
dayfunc<-function(value){
  if(value==1){
    weekdayfinal<-"Sunday"
  }
  if(value==2){
    weekdayfinal<-"Monday"
  }
  if(value==3){
    weekdayfinal<-"Tuesday"
  }
  if(value==4){
    weekdayfinal<-"Wednesday"
  }
  if(value==5){
    weekdayfinal<-"Thursday"
  }
  if(value==6){
    weekdayfinal<-"Friday"
  }
  if(value==7){
    weekdayfinal<-"Saturday"
  }
  day<-weekdayfinal
}
predictfunc<-function(dataframe){
  predictionlocation1<-dataframe
  timepredict<-predict.lm(week1elapsedtime, predictionlocation1, se.fit = TRUE)
  totelapsed<-exp(timepredict$fit)
  timeprediction<-as.data.frame(totelapsed)
  longpredict<-predict.lm(week1long, timeprediction, se.fit = TRUE)
  finallong<-data.frame(long=exp(longpredict$fit[1]), totelapsed=totelapsed)
  latpredict<-predict.lm(week1lat, finallong, se.fit = TRUE)
  finallat<-as.data.frame(exp(latpredict$fit))
  finalloc<-data.frame(long=finallong[,1], lat=finallat[,1], totelapsed=totelapsed)
  predictionvalue<-dataframe$time_stamp[1]+totelapsed
  dayvalue<-wday(predictionvalue)
  day<-dayfunc(dayvalue)
  predictionvalue<-data.frame(day=day, 
                              time_stamp=predictionvalue[1], long=finalloc$long, lat=finalloc$lat)
  if (totelapsed<913 |totelapsed>1697){
    print("Do Not Bomb This Trip")
  }
  else {
    return(predictionvalue)
  }
}

#Read in CSV file and plug in resulting list to function below
finalfunc<-function(csvdata){
  initialdat<-data.frame(long=NA, lat=NA, time_stamp=NA)
  initialdata<-data.frame()
  i<-1
  for (i in 1:length(csvdata[["features"]])){
    initialdat$long[1]<-csvdata[["features"]][[i]][["geometry"]][["coordinates"]][[1]]
    initialdat$lat[1]<-csvdata[["features"]][[i]][["geometry"]][["coordinates"]][[2]]
    initialdat$time_stamp[1]<-csvdata[["features"]][[i]][["properties"]][["time"]]
    initialdata<-rbind(initialdata, initialdat)
  }
  initialdata$time_stamp<-as.POSIXct(initialdata$time_stamp, format="%Y-%m-%dT%H:%M:%S")
  spat_df <- SpatialPointsDataFrame(coords=initialdata[, c("long", "lat")],
                                    data=initialdata['time_stamp'], 
                                    proj4string=CRS("+proj=longlat +datum=WGS84"))
  utm_df <- spTransform(spat_df, CRSobj = "+proj=utm +zone=12 +datum=WGS84")
  utm_coords <- coordinates(utm_df)
  utmdataframe<-c()
  utmdataframe$long<-utm_coords[,1]
  utmdataframe$lat<-utm_coords[,2]
  utmdataframe<-as.data.frame(utmdataframe)
  utmdataframe$time_stamp<-utm_df@data[["time_stamp"]]
  return(predictfunc(utmdataframe))
}

finalfunc(testgps)
### end
##### Temp ####
temp<-read.csv("~/Desktop/2338340.csv")

temp<-temp[temp$STATION=="USR0000MMIS",]
temp<-temp[!temp$DATE=="2020-08-17",]
temp<-temp[!temp$DATE=="2020-08-23",]
temp<-temp[!temp$DATE=="2020-08-29",]
temp<-temp[!temp$DATE=="2020-08-30",]
### end

##### Speed ####
df2<-utmdatafull
df2[which(df2$elapsedtime > 120),]<-NA
dt<-mean(df2$elapsedtime, na.rm=T)
dt

i<-1
for (i in 1:nrow(utmdatafull)){
  if(utmdatafull$totelapsed[i]==0){
    utmdatafull$velx[i]<-0
    utmdatafull$vely[i]<-0
  }
  if(utmdatafull$elapsedtime[i]==0){
    utmdatafull$velx[i]<-0
    utmdatafull$vely[i]<-0
  }
  else{
  utmdatafull$velx[i]<-(utmdatafull$long[i]-utmdatafull$long[i-1])/utmdatafull$elapsedtime[i]
  utmdatafull$vely[i]<-(utmdatafull$lat[i]-utmdatafull$lat[i-1])/utmdatafull$elapsedtime[i]
  }}


i<-1
for (i in 1:nrow(utmdatafull)){
  utmdatafull$speed[i]<-sqrt((utmdatafull$velx[i])^2+(utmdatafull$vely[i])^2)
  
}

speedfunc<-function (x){
  speed<-sqrt(x[[7]]^2+x[[8]]^2)
}
speed<-speedfunc(utmdatafull)
speed
speed<-apply(utmdatafull,1,speedfunc)


Fmat<-matrix(c(1, 0, dt, 0,
               0, 1, 0, dt), byrow=TRUE, ncol=4)

gps_variance <- 20^2
v_mat <- matrix(c(gps_variance, 0,
                  0, gps_variance), byrow=TRUE, ncol=2)
dt<-mean(df2$elapsedtime, na.rm=T) #8.15
g_mat <- matrix(c(1, 0, dt, 0,
                  0, 1, 0, dt,
                  0, 0, 1, 0,
                  0, 0, 0, 1), byrow=TRUE, ncol=4)
avg_walk_speed_m_per_sec <- 1.4  # https://en.wikipedia.org/wiki/Walking
dlm_spec <- dlm(
  FF= Fmat,
  GG= g_mat,
  V = v_mat,
  W = diag(c(5, 5, 1, 1)^2),
  m0 = matrix(c(utm_coords[1, ], rep(avg_walk_speed_m_per_sec / dt, 2)),
              ncol=1), # A vector by R defaults is a k by 1 matrix
  C0 = diag(rep(10^2, 4)))

dlm_filter_mod <- dlmFilter(utm_coords, dlm_spec)
dlm_smooth_mod <- dlmSmooth(dlm_filter_mod)

speed<-c()

speeddf<-dlm_smooth_mod$s[-1,]

speedfunc<-function (x){
  speed<-sqrt(x[[3]]^2+x[[4]]^2)
}

speed<-apply(speeddf,1,speedfunc)
### end


##### Combine speed/temp ####
utmdatafull$speed<-speed
df1<-utmdatafull
y<-unique(mday(utmdatafull$time_stamp))
date<-c()
speedvec<-c()
i<-1
for (i in 1:length(y)){
  date1<-y[i]
  df3<-df1$speed[which(mday(df1$time_stamp)== y[i])]
  speed1<-mean(df3)
  date<-rbind(date, date1)
  speedvec<-rbind(speedvec, speed1)
}
tempspeed<-c()
nrow(temp)
tempspeed<-cbind(date, speedvec, temp$TAVG)
tempspeed<-data.frame(Date=date, Speed=speedvec, TAvg=temp$TAVG)

ggplot(tempspeed, aes(TAvg, Speed))+
  geom_point()+
  theme_minimal()+
  xlab("Average Temperature (C)")+
  ylab("Speed (m/s)")+
  ggtitle("Average Daily Speed and Temperature")
  geom_smooth(method="lm")
cor(tempspeed$TAvg, tempspeed$Speed)
### end


##### Plots ####
quantile(utmtrain$totelapsed)
meansmall<-mean(utmtrain$totelapsed[!utmtrain$totelapsed<626 & !utmtrain$totelapsed>2000])
sd<-sd(utmtrain$totelapsed[utmtrain$totelapsed>626&utmtrain$totelapsed<2000])


ggplot(utmdatafull[utmdatafull$totelapsed==0,], aes(wday(time_stamp), hour(time_stamp), color=as.factor(session)))+
  geom_point()+
  xlab("Weekday")+
  ylab("Hour of Day")+
  ggtitle("Hour and Day of Departure")+
  theme_minimal()+
  scale_color_discrete(name="Session Number")

ggplot(utmdatafull, aes(long, lat))+
  geom_point()+
  theme_minimal()+
  ggtitle("Location Data for Entire Period")

ggplot(utmdatafull[mday(utmdatafull$time_stamp)>23,], aes(time_stamp, lat, color=as.factor(session)))+
  geom_point()+
  theme_minimal()+
  xlab("Time Stamp")+
  ylab("Latitude")+
  ylim(5193950, 5197050)+
  ggtitle("Latitude vs. Time for Entire Data")+
  scale_color_discrete(name="Session Number")

ggplot(utmdatafull[mday(utmdatafull$time_stamp)>23,], aes(time_stamp, long, color=as.factor(session)))+
  geom_point()+
  theme_minimal()+
  xlab("Time Stamp")+
  ylab("Longitude")+
  ylim(271400, 272550)+
  ggtitle("Longitude vs. Time for Entire Data")+
  scale_color_discrete(name="Session Number")

ggplot(utmtrain, aes(long, lat))+
  geom_point()+
  theme_minimal()+
  ggtitle("Location Data for Cleaned Data")

ggplot(utmtrain, aes(totelapsed, long))+
  geom_point()+
  theme_minimal()+
  xlim(626,2000)+
  xlab("Elapsed Time")+
  ylab("Longitude")+
  ggtitle("Longitude and Elapsed Time for Cleaned Data")+
  geom_vline(xintercept=meansmall, color="red")+
  geom_vline(xintercept=(meansmall+sd), color="green")+
  geom_vline(xintercept=(meansmall-sd), color="green")

ggplot(utmtrain, aes(totelapsed, lat))+
  geom_point()+
  theme_minimal()+
  xlim(626,2000)+
  xlab("Elapsed Time")+
  ylab("Latitude")+
  ggtitle("Latitude and Elapsed Time for Cleaned Data")+
  geom_vline(xintercept=meansmall, color="red")+
  geom_vline(xintercept=(meansmall+sd), color="green")+
  geom_vline(xintercept=(meansmall-sd), color="green")


ggplot(utmtrain, aes(lat, long, color=totelapsed))+
  geom_point()+
  xlab("Latitude")+
  ylab("Longitude")+
  theme_minimal()+
  ggtitle("Latitude/Longitude and Elapsed Time")+
  scale_color_continuous(name="Elapsed Time")

ggplot(utmtrain, aes(totelapsed,long))+
  geom_point()+
  theme_minimal()+
  xlab("Elapsed Time")+
  ylab("Longitude")+
  ggtitle("Longitude and Elapsed Time")

ggplot(utmtrain, aes(totelapsed,lat))+
  geom_point()+
  theme_minimal()+
  xlab("Elapsed Time")+
  ylab("Longitude")+
  ggtitle("Longitude and Elapsed Time")

p1<-ggplot(utmtrain[utmtrain$session==10,], aes(long,lat))+
  geom_point()+
  theme_minimal()+
  xlab("Longitude")+
  ylab("Latitude")+
  geom_point(data=finalloc,color="red")+
  geom_point(color="green", aes(x=271781, y= 5195721))+
  geom_point(color="green", aes(x=271781, y= 5195260))+
  geom_point(color="green", aes(x=272059.9, y= 5195491))+
  geom_point(color="green", aes(x=271502.4, y= 5195491))

p2<-ggplot(utmtrain[utmtrain$session==10,], aes(totelapsed,long))+
  geom_point()+
  theme_minimal()+
  ylab("Longitude")+
  xlab("Elapsed Time")+
  geom_point(data=finalloc,color="red")+  
  geom_point(color="green", aes(x=1307.548, y= 272059.9))+
  geom_point(color="green", aes(x=1307.548, y= 271502.4))+
  geom_point(color="green", aes(x=1447.362, y= 271781))+
  geom_point(color="green", aes(x=1181.24, y= 271781))


p3<-ggplot(utmtrain[utmtrain$session==10,], aes(totelapsed,lat))+
  geom_point()+
  theme_minimal()+
  ylab("Latitude")+
  xlab("Elapsed Time")+
  geom_point(data=finalloc,color="red")+
  geom_point(color="green", aes(x=1307.548, y= 5195721))+
  geom_point(color="green", aes(x=1307.548, y= 5195260))+
  geom_point(color="green", aes(x=1447.362, y= 5195491))+
  geom_point(color="green", aes(x=1181.24, y= 5195491))

ggarrange(p1,p2,p3, nrow=2, ncol=2)
### end
