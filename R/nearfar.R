
#Function to compute rank-based Mahalanobis distance matrix between each pair based on observed covariates (make "near"); Paul Rosenbaum is author
#X is a matrix of observed covariates with n rows (observations) and p columns (variables)
smahal = function(X){
  X<-as.matrix(X)
  n<-dim(X)[1]
  k<-dim(X)[2]
  for (j in 1:k) X[,j]<-rank(X[,j]) #compute on ranks
  cv<-cov(X)
  vuntied<-var(1:n)
  rat<-sqrt(vuntied/diag(cv))
  cv<-diag(rat)%*%cv%*%diag(rat)
  out<-matrix(NA,n,n)
  icov<-ginv(cv)
  for (i in 1:n) out[i,]<-mahalanobis(X,X[i,],icov,inverted=T)
  out
}

#Prioritize importance of matches on certain observed covariates
#Distmat is an object of class distance matrix
#variable is named variable from the measured covariates to update the distance matrix to prioritize "nearness" of that variable
#tolerance is how much priority to attach to the named variable; values close to 0 indicate HIGH priority
calipers = function(distmat,variable,tolerance=0.2){
  sd.distmat = median(sd(distmat))
  sd.var = sd(variable)
  penalty = matrix(variable,nrow=length(variable),ncol=length(variable))
  penalty = penalty - t(penalty) #how far apart each pair of individuals are with regard to variable of interest
  distmat[abs(penalty)> tolerance*sd.var]  <- distmat[abs(penalty)> tolerance*sd.var]  + abs((penalty/sd.var)[abs(penalty)> tolerance*sd.var]*.5*sd.distmat) #add pentalty of abs((distance/sd.var)*.5*sd.distmat) to parts of distance matrix that are far enough apart with regard to variable of interest
  distmat
}

#Function to carry out nearfar match
#X is a vector of measured covariates (with column names) on which to make "near" in the matching
#imp.var is a list of (up to 5) named variables to prioritize in the "near" matching
#tol.var is list of (up to 5) tolerances attached to the prioritized variables where 0 is highest priority
#sinks is what percent of the data to match to sinks (and thus remove) if desired; default is 0
#cutpoint is value below which individuals are too similar on IV; increase to make individuals more "far" in match
#IV is a vector of instrumental variables on which to make "far" in the matching
#Default settings yield a "near" match on only observed covariates in X

matches = function(X,imp.var=NA,tol.var=NA,sinks=0,IV=NA,cutpoint=NA){
  cols.w.var <-  which(apply(X,2,var)>0)#Only give the columns with variation
  distmat=smahal(X[,cols.w.var]);#Create distance matrix from observed covariates WITH VARIATION
  distmat=round(distmat*1000);#Make into large integer distance so sinks are smallest.

#Prioritizing the observed covariate match
  if (sum(is.na(imp.var)==0)>0) {
  if (length(imp.var)>0) {
      j1 = calipers(distmat,X[,names(X)==imp.var[1]],tolerance=tol.var[1])
      distmat = j1
  }
  if (length(imp.var)>1) {
      j2 = calipers(j1,X[,names(X)==imp.var[2]],tolerance=tol.var[2])
      distmat = j2
  }
  if (length(imp.var)>2) {
      j3 = calipers(j2,X[,names(X)==imp.var[3]],tolerance=tol.var[3])
      distmat = j3
  }
  if (length(imp.var)>3) {
      j4 = calipers(j3,X[,names(X)==imp.var[4]],tolerance=tol.var[4])
      distmat = j4
  }
  if (length(imp.var)>4) {
      j5 = calipers(j4,X[,names(X)==imp.var[5]],tolerance=tol.var[5])
      distmat = j5
  }
}

#Update distance matrix to make matches "far" on IV
  if (sinks>0) {
  distmat.adjust <- min(sd(distmat)) #important for penalty later
  mat.IV <- matrix(IV,nrow=length(IV),ncol=length(IV)) #continuous IV
  dist.IV <- abs(mat.IV - t(mat.IV)) #pairwise IV discrepancies
  #cutpoint is the value of IV below which individuals are encouraged
  distmat[dist.IV<cutpoint]  <-   distmat[dist.IV<cutpoint]  + (max(dist.IV[dist.IV<cutpoint])-dist.IV[dist.IV<cutpoint])*distmat.adjust #Penalize individuals at similar level of IV by increasing distance between them in distance matrix

  #Create sinks if desired
  size <- dim(distmat)[2]
  #sinks=0.1
  num.sinks <- size*sinks  #The "sinks" variable tells you what fraction of the individuals should be matched to sinks, i.e., removed from analysis
  num.sinks <- 2*ceiling((size+num.sinks)/2) - size        #Make sure we have an even number of (individuals + sinks) to match
  total <- size + num.sinks
  distmat <- cbind(rbind(distmat,matrix(0,nrow=num.sinks,ncol=size)),matrix(0,nrow=total,ncol=num.sinks));     #Add all of the sinks; they are are 0 distance from every other observation, but there are only so many of them, so not everyone can be matched to a sink; only the most "consistently discordant" individuals are matched to sinks
  if(num.sinks>0){distmat[(size+1):total,(size+1):total] <- max(distmat)*3} #If sinks, make sure sinks don't match to sinks
}
  #Perform the matching on distance matrix
  distmat=round(distmat);
  distmat <- distancematrix(distmat) #Make into distancematrix object
  matching=nonbimatch(distmat); #match on observed covs; optionally on IV
  encouraged.numbers = matching$halves[,2] #define initial vectors
  discouraged.match = matching$halves[,4]
  #Redefine vectors if strengthening IV
  if (sinks>0) {
  IV.orig=c(IV,rep(min(IV)-1,num.sinks));     #Be sure to add SMALL IV values for the sinks
  IV.matched=IV[matching$matches[,4]] #NAs for the sinks
  encouraged=IV.orig<IV.matched;  #Figuring out which one in the pair is encouraged; orig for sinks is 1; matched for sinks is NA - there are going to be NAs because of the sinks
  if (num.sinks > 0) encouraged[(size+1):total]=FALSE;
  encouraged[is.na(IV.matched)]=FALSE; #The NAs tell you that there's been an individual matched to a sink.

  #Some IVs will be the same within a pair; let half be encouraged, half be discouraged and if odd number just drop last pair
  same.IV <- sum(IV.orig==IV.matched,na.rm=TRUE)
  IV.matched[is.na(IV.matched)] <- IV.orig[is.na(IV.matched)]+1
  if(same.IV>0){
    encouraged[which(IV.orig==IV.matched)] <- sample(c(rep(TRUE,floor(same.IV/2)),rep(FALSE,ceiling(same.IV/2))),size=same.IV)
}
  encouraged.numbers=which(encouraged); #numbers of encouraged
  discouraged.match=matching$matches[encouraged.numbers,4];
}
  return(cbind(encouraged.numbers,discouraged.match))
}

#Compute absolute standardized difference for one variable
#x is a vector of values of a continous or binary variable
#match is a two column matrix where the first column is the index of an "encouraged" individual and the second column is the index of the corresponding "discouraged" individual from the pair matching
asd = function(x,match) {
    enc = x[match[,1]]
    unenc = x[match[,2]]
    m.enc = mean(enc,na.rm=TRUE)
    m.unenc = mean(unenc,na.rm=TRUE)
    sd2 = sqrt(var(enc,na.rm=TRUE)/2+var(unenc,na.rm=TRUE)/2) #change this
    x2 = x[is.na(x)==0]
    t = as.numeric(names(table(x2)))
    if (length(t)==2 & t[1]==0 & t[2]==1) {
        sd2 = sqrt(m.enc*(1-m.enc)/2+m.unenc*(1-m.unenc)/2) #change for binary variable
    }
    asd = abs(m.enc-m.unenc)/sd2
    output.all = c(m.enc,m.unenc,asd)
    return(output.all)
}

#Function to evaluate nearfar match
#dta is a data frame wherein the rows are individuals and the variables are outcome, exposure, IV, and covariates
#match is a two column matrix where the first column is the index of an "encouraged" individual and the second column is the index of the corresponding "discouraged" individual from the pair matching
summ.matches = function(dta,match) {
    npairs = dim(match)[1]
res = t(apply(dta,2,function(x) asd(as.numeric(x),match)))
colnames(res) = c(paste("Encouraged n=",npairs,sep=""),paste("Discouraged n=",npairs,sep=""),"Abs St Dif")
#print(paste(npairs,"pairs of observations",sep=" "))
return(res)
}

#Function to optimize nearfar match
#dta is data frame wherein first column is outcome, second column is exposure, third column is IV, and fourth through last columns are measured confounders
#trt.bin is an indicator of a binary exposure variable; default is FALSE; if true, the sink and cutpoint parameters in the near-far match are optimized based on residual deviance for the IV
#imp.var is a list of (up to 5) named variables to prioritize in the "near" matching
#tol.var is list of (up to 5) tolerances attached to the prioritized variables where 0 is highest priority
#adjust.IV is whether or not to include measured confounders in exposure~IV model that is optimized
#sink.range is a two element vector of (min, max) for range of sinks over which to optimize in the near-far match; default (0,0.5) such that maximally 50% of observations can be removed
#cutp.range is a two element vector of (min, max) for range of cutpoints (how far apart the IV will become) over which to optimize in the near-far match; default is min=one SD of IV, max=range of IV
#max.time.seconds is how long to let the optimization algorithm run; default is 300 seconds = 5 minutes
opt.nearfar = function(dta,trt.bin=FALSE,imp.var=NA,tol.var=NA,adjust.IV=TRUE,sink.range=c(0,0.5),cutp.range=NA,max.time.seconds=300) {
    df = dta[,-1] #remove outcome
    #expo = df[,1]
    IV = df[,2]
    IV2 = as.numeric(IV)
    if (is.na(sum(cutp.range))==1) {cutp.range = c(sd(IV2),max(IV2)-min(IV2))} #make range of cutpoints [one SD of IV, range of IV]
    X = df[,3:dim(df)[2]]
    if (length(dim(X))==0) {X2 = matrix(as.numeric(as.matrix(X)),length(X),1)} #for 1 measured covariate
    else if (length(dim(X))==2) {X2 = matrix(as.numeric(as.matrix(X)),dim(X)[1],dim(X)[2])} #for >1 measured covariate
#Function to be minimized
if (trt.bin==FALSE & adjust.IV==FALSE) {
    findf = function(x) {
    pars = x
    r = matches(X=X2,imp.var,tol.var,sinks=pars[1],IV=IV2,cutpoint=pars[2])
#Does near-far match improve F-statistic?
    df$match.NF = NA #add last column to df
    df$match.NF[r[,1]] = 1 #encouraged
    df$match.NF[r[,2]] = 0 #discouraged match
    df2 = df[is.na(df$match.NF)==0,] #6000 or (1-sinks)*N sample size
    df2 = df2[,-dim(df2)[2]]
    names(df2)[1] = "expo"
    names(df2)[2] = "IV"
    m2 = lm(expo~IV,data=df2)
    j = anova(m2)
    f = j[which(rownames(j)=="IV"),which(colnames(j)=="F value")]
    -f
 }
}
else if (trt.bin==FALSE & adjust.IV==TRUE) {
    findf = function(x) {
    pars = x
    r = matches(X=X2,imp.var,tol.var,sinks=pars[1],IV=IV2,cutpoint=pars[2])
#Does near-far match improve F-statistic?
    df$match.NF = NA #add last column to df
    df$match.NF[r[,1]] = 1 #encouraged
    df$match.NF[r[,2]] = 0 #discouraged match
    df2 = df[is.na(df$match.NF)==0,] #6000 or (1-sinks)*N sample size
    df2 = df2[,-dim(df2)[2]]
    names(df2)[1] = "expo"
    names(df2)[2] = "IV"
    df3 = df2[,c(1,3:dim(df2)[2],2)] #re-order for partial F
    m2 = lm(expo~.,data=df3)
    j = anova(m2)
    f = j[which(rownames(j)=="IV"),which(colnames(j)=="F value")]
    -f
 }
}
else if (trt.bin==TRUE & adjust.IV==FALSE) {
    findf = function(x) {
    pars = x
    r = matches(X=X2,imp.var,tol.var,sinks=pars[1],IV=IV2,cutpoint=pars[2])
#Does near-far match improve F-statistic?
    df$match.NF = NA #add last column to df
    df$match.NF[r[,1]] = 1 #encouraged
    df$match.NF[r[,2]] = 0 #discouraged match
    df2 = df[is.na(df$match.NF)==0,] #6000 or (1-sinks)*N sample size
    df2 = df2[,-dim(df2)[2]]
    names(df2)[1] = "expo"
    names(df2)[2] = "IV"
    m3 = glm(expo~IV,data=df2,family=binomial)
    j = anova(m3,test="Chisq")
    d = j[which(rownames(j)=="IV"),which(colnames(j)=="Deviance")]
    -d
 }
}
else if (trt.bin==TRUE & adjust.IV==TRUE) {
    findf = function(x) {
    pars = x
    r = matches(X=X2,imp.var,tol.var,sinks=pars[1],IV=IV2,cutpoint=pars[2])
#Does near-far match improve F-statistic?
    df$match.NF = NA #add last column to df
    df$match.NF[r[,1]] = 1 #encouraged
    df$match.NF[r[,2]] = 0 #discouraged match
    df2 = df[is.na(df$match.NF)==0,] #6000 or (1-sinks)*N sample size
    df2 = df2[,-dim(df2)[2]]
    names(df2)[1] = "expo"
    names(df2)[2] = "IV"
    df3 = df2[,c(1,3:dim(df2)[2],2)] #re-order for partial d
    m3 = glm(expo~.,data=df3,family=binomial)
    j = anova(m3,test="Chisq")
    d = j[which(rownames(j)=="IV"),which(colnames(j)=="Deviance")]
    -d
 }
}
#Optimize using GenSA (simulated annealing)
    out <- GenSA(lower=c(sink.range[1],cutp.range[1]),upper=c(sink.range[2],cutp.range[2]),fn=findf,control=list(max.time=max.time.seconds))
    val = out$value
    pct.sink = out$par[1]
    cutp = out$par[2]
    n.calls = out$counts
    t = matches(X=X2,imp.var,tol.var,sinks=pct.sink,IV=IV2,cutpoint=cutp)
    summ = summ.matches(dta=dta,match=t)
    output.all = list(n.calls=n.calls,sink.range=sink.range,cutp.range=cutp.range,pct.sink=pct.sink,cutp=cutp,maxF=-val,match=t,summ=summ)
    return(output.all) #get %sinks, cutpoint, and matches
}

#Compute test statistic as a function of lambda
stat.lambda = function(lambda,dta,match) {
    R = data.frame(R.t=dta[match[,1],1],R.c=dta[match[,2],1]) #outcome
    D = data.frame(D.t=dta[match[,1],2],D.c=dta[match[,2],2]) #treatment
    V = (R$R.t-lambda*D$D.t)-(R$R.c-lambda*D$D.c)
    T0 = mean(V)
    I = length(V)
    S0 = sqrt(sum((V-T0)^2)/(I*(I-1)))
    return(T0/S0)
}

#Function for inference on the effect ratio
eff.ratio = function(dta,match,alpha) {
    R = data.frame(R.t=dta[match[,1],1],R.c=dta[match[,2],1]) #outcome
    D = data.frame(D.t=dta[match[,1],2],D.c=dta[match[,2],2]) #treatment
    est1 = mean(R$R.t-R$R.c)/mean(D$D.t-D$D.c)
    est2 = tryCatch(uniroot(function(x) stat.lambda(x,dta,match),lower=-10000,upper=10000,tol=1e-9)$root,error=function(e){})
    v1 = NA
    v1 = tryCatch(uniroot(function(x) stat.lambda(x,dta,match)-qnorm(1-alpha/2),lower=-10000,upper=10000,tol=1e-9)$root,error=function(e){})
    v2 = NA
    v2 = tryCatch(uniroot(function(x) stat.lambda(x,dta,match)+qnorm(1-alpha/2),lower=-10000,upper=10000,tol=1e-9)$root,error=function(e){})
    lower = min(v1,v2)
    upper = max(v1,v2)
    output.all = list(est.emp=est1,est.HL=est2,lower=lower,upper=upper)
    return(output.all)
}

#eff.ratio(dta,match)

#Example for paper
#l = seq(-10,10,0.01)
#stat = as.numeric(lapply(l,function(x) stat.lambda(x,dta,match)))
#pdf("/Users/jrigdon/Box Sync/Rigdon/Sanjay/Paper1/nf_2016-08-09.pdf")
#plot(l,stat,xlab="Effect Ratio",ylab="Statistic")
#dev.off()

#package.skeleton(name="nearfar", code_files="nearfar.R")




