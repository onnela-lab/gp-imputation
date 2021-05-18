library(circular)

####################################################################################
## part 1: mother traj
n = 2160
## p(p|p) = 0.85, p(f|p) = 0.15
## p(f|f) = 0.7, p(p|f) = 0.3
## sp = runif(1,1.4)
## t = runif(20,60)
## dir_0 = rvonmises(n, 0, kappa,control.circular=list(units="degrees"))
## dir_i = rvonmises(n, i-1, kappa,control.circular=list(units="degrees"))

s = 1
dir = circular(0)
l_c = c(0,0)
trace = l_c
t_seq = c(0)
s_seq = c()
a_seq = c()
sp_seq = c()
t_diff = c()
for(i in 1:n){
  if(s==1){
    new_s = rbinom(1,1,0.7)
  }
  if(s==0){
    new_s = rbinom(1,1,0.15)
  }
  s = new_s
  s_seq = c(s_seq,s)
  if(s==1){
    new_dir = 0.5*(rvonmises(1, circular(dir),10,control.circular=list(units="radians"))+dir)
    dir = new_dir
    a_seq = c(a_seq,as.numeric(dir))
    sp = runif(1,1,1.4)
    sp_seq = c(sp_seq,sp)
    t = runif(1,20,60)
    mov = c(sp*t*cos(as.numeric(dir)),sp*t*sin(as.numeric(dir)))
    l_c = l_c + mov
    trace = rbind(trace,l_c)
    t_seq = c(t_seq,t_seq[length(t_seq)]+t)
  }
  if(s==0){
    t = runif(1,20,60)
    trace = rbind(trace,l_c)
    t_seq = c(t_seq,t_seq[length(t_seq)]+t)
  }
  t_diff = c(t_diff,t)
}

trace0 = trace
t_seq0 = t_seq
plot(trace0[,1],trace0[,2],"l")

###################################################################
## part2: offspring traj
n = 2160
k = function(t,x,y,t_seq,trace){
  k1 = exp(-abs(t-t_seq)/(60*10))
  d = sqrt((x-trace[,1])^2+(y-trace[,2])^2)
  k2 = exp(-d/200)
  return(0.5*k1+0.5*k2)
}

s = 1
dir = circular(0)
l_c = c(0,0)
trace = matrix(l_c,nrow=1)
t_seq = c(0)
for(i in 1:n){
  x = trace[i,1]
  y = trace[i,2]
  t = t_seq[i]
  k0 = k(t,x,y,t_seq0[1:n],trace0[1:n,])
  if(s==1){
    p = 0.5*(0.7 + sum(k0[s_seq==1])/sum(k0))
    new_s = rbinom(1,1,p)
  }
  if(s==0){
    p = 0.5*(0.15 +sum(k0[s_seq==1])/sum(k0)) 
    new_s = rbinom(1,1,p)
  }
  s = new_s
  if(s==1){
    guess_dir = circular(sum(k0[s_seq==1]*a_seq)/sum(k0[s_seq==1]))
    new_dir = as.numeric(rvonmises(1, 0.5*(circular(dir)+guess_dir),10,control.circular=list(units="radians")))
    dir = new_dir
    guess_sp = sum(k0[s_seq==1]*sp_seq)/sum(k0[s_seq==1])
    sp = rnorm(1,guess_sp,0.02)
    guess_t = sum(k0[s_seq==1]*t_diff[s_seq==1])/sum(k0[s_seq==1])
    t = rnorm(1,guess_t,2)
    mov = c(sp*t*cos(dir),sp*t*sin(dir))
    l_c = l_c + mov
    trace = rbind(trace,l_c)
    t_seq = c(t_seq,t_seq[length(t_seq)]+t)
  }
  if(s==0){
    guess_t = sum(k0[s_seq==0]*t_diff[s_seq==0])/sum(k0[s_seq==0])
    t = rnorm(1,guess_t,2)
    trace = rbind(trace,l_c)
    t_seq = c(t_seq,t_seq[length(t_seq)]+t)
  }
}

trace1 = trace
t_seq1 = t_seq
plot(trace1[,1],trace1[,2],"l")


########################################################################################
## construct full data
full_data = c()
t_seq1 = round(t_seq1,0)
if(max(t_seq1)>86400){
  index = t_seq1<86400
  t_seq1 = t_seq1[index]
  trace1 = trace1[index,]
}

if(max(t_seq1)<86400){
  t_seq1 = c(t_seq1,86400)
  trace1 = rbind(trace1,trace1[nrow(trace1),])
}

j = 1
for(i in 0:86400){
  print(i/86400)
  if(sum(t_seq1==i)==1){
    j = (1:nrow(trace1))[t_seq1==i]
    full_data = rbind(full_data,c(i,trace1[j,]))
  }
  else{
    t0 = t_seq1[j]
    t1 = t_seq1[j+1]
    lat = (t1-i)/(t1-t0)*trace1[j,1]+(i-t0)/(t1-t0)*trace1[j+1,1]
    lon = (t1-i)/(t1-t0)*trace1[j,2]+(i-t0)/(t1-t0)*trace1[j+1,2]
    full_data = rbind(full_data,c(i,lat,lon))
  }
}

full_data = data.frame(matrix(as.numeric(full_data),ncol=3))
full_data$accuracy = 20
colnames(full_data) = c("timestamp","latitude","longitude","accuracy")

###############################################################################
## construct observed data
obsdata = c()
start_t = 0
t_seq = 0:86400
n = 86400/(10*60)
for(i in 1:n){
  end_t = start_t + 600
  obsdata = rbind(obsdata,full_data[t_seq>=start_t&t_seq<=end_t,])
  start_t = start_t + rnorm(1,120,1.5)*60
}
plot(obsdata$latitude,obsdata$longitude)

write.csv(obsdata,"day3_obs.csv",row.names = F)
write.csv(full_data,"day3_full.csv",row.names = F)

############################################################################


data1 = read.csv("day1_obs.csv")
data2 = read.csv("day2_obs.csv")
data3 = read.csv("day3_obs.csv")
data1 = data1[data1$timestamp<3600|data1$timestamp>25200,]
data2 = data2[data2$timestamp<32400|data2$timestamp>54000,]
data3 = data3[data3$timestamp<61200|data3$timestamp>82800,]
data2$timestamp = data2$timestamp + 86400
data3$timestamp = data3$timestamp + 86400*2
obs = rbind(data1,data2,data3)

## convert to lat/lon
obs$latitude = obs$latitude/11119.5*0.1+42
obs$longitude = obs$longitude/8263.3*0.1-71
write.csv(obs,"vonmises_obs.csv",row.names = F)
