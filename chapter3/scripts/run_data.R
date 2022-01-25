## Usage: "Rscript run_data.R <inputfile> <outputroot> <Fst> <seed>" for integer seeds to run the stan models
## inputfile should be a "snpRes" file from BayesS
## outputroot is used for output
## Fst is set for distance from original population. 0.15 is a good value
## seed is for repeatability. seed=1 is ok

args=commandArgs(TRUE)
if(length(args)!=4){
    stop("Usage: Rscript run_data.R <inputfile> <name> <Fst> <seed>")
}

## Parameters we might like to set; you have to change these by hand!
iter=10000
thin=50
thresh=0.001
ncores=2

inputfile=args[[1]]
outputroot=args[[2]]
Fst=as.numeric(args[[3]])
seed=as.numeric(args[[4]])
set.seed(seed)

print(paste("Using Data:", inputfile))
print(paste("Using outputroot:", outputroot))
print(paste("Using Fst:", Fst))
print(paste("Using seed:", seed))

cpgtab=read.table(inputfile,header=T, row.names=1,as.is=T)
cpgtab=cpgtab[cpgtab$PIP>thresh,] # From EDA of the single CPG I saw...
cpgtab=cpgtab[order(cpgtab$PIP,decreasing=T),]

## Functions
data_as_stan_list=function(x,Fst=0.15){
    list(J=dim(x)[1],
         beta=x$Effect,
         f=x$Freq,
         w=x$PIP,
         s=x$SE,
         Fst=Fst)
}
getS=function(stanres){
    svals=as.numeric(extract(stanres,"S")[[1]])
    c(mean=mean(svals),
      sd=sd(svals),
      quantile(svals,c(0.05,0.25,0.5,0.75,0.95)))
}
## set Fst now
cpgdat=data_as_stan_list(cpgtab,Fst)

########
## STAN
library("rstan")
options(mc.cores=ncores)

sm0=readRDS("model0.RDS")
sm5=readRDS("model5.RDS")
sm6=readRDS("model6.RDS")
sm7=readRDS("model7.RDS")
res0<-sampling(sm0,data=cpgdat,chains=ncores,iter=iter,thin=thin)
res5<-sampling(sm5,data=cpgdat,chains=ncores,iter=iter,thin=thin)
res6<-sampling(sm6,data=cpgdat,chains=ncores,iter=iter,thin=thin)
res7<-sampling(sm7,data=cpgdat,chains=32,iter=iter,thin=thin)

allres=list(model0=res0,
            model5=res5,
            model6=res6,
            model7=res7)
resdf=t(sapply(allres,getS))

write.table(resdf,file=paste0(outputroot,".csv"),quote=FALSE,sep=",")
saveRDS(allres,file=paste0(outputroot,".RDS"))
