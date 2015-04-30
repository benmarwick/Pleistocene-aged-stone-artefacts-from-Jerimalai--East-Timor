graphics.off()
# rm(list=ls(all=TRUE))
source("code/openGraphSaveGraph.R")
source("code/plotPost.R")
fileNameRoot="ANOVAonewayJagsSTZ" # for constructing output filenames
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
                       # A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.

modelstring = "
model {
  for ( i in 1:Ntotal ) {
    y[i] ~ dnorm( mu[i] , tau )
    mu[i] <- a0 + a[x[i]]
  }
  #
  tau <- 1 / pow( sigma , 2 )
  sigma ~ dunif(0,10) # y values are assumed to be standardized
  #
  a0 ~ dnorm(0,0.001) # y values are assumed to be standardized
  #
  for ( j in 1:NxLvl ) { a[j] ~ dnorm( 0.0 , atau ) }
  atau <- 1 / pow( aSD , 2 )
  aSD ~ dgamma(1.01005,0.1005) # mode=0.1,sd=10.0
  # Convert a0,a[] to sum-to-zero b0,b[] :
  for ( j in 1:NxLvl ) { m[j] <- a0 + a[j] } 
  b0 <- mean( m[1:NxLvl] )
  for ( j in 1:NxLvl ) { b[j] <- m[j] - b0 }
}
" # close quote for modelstring
writeLines(modelstring,con="model.txt")

#--------
# THE DATA.

# Specify data source:
dataSource = c( "data" )[1]
# Load the data:

if ( dataSource == "data" ) {
  fileNameRoot = paste( fileNameRoot , dataSource , sep="" )
  datarecord = data 
  y = as.numeric(datarecord[,2]) # numeric measurements
  Ntotal = length(datarecord[,2]) # numeric groups
  x = as.numeric(datarecord[,1]) # numeric groups
  xnames = levels(as.factor(datarecord[,1]))
  NxLvl = length(unique(datarecord[,1]))
  normalize = function( v ){ return( v / sum(v) ) }
  contrastList = list( 
    # note that these are custom interactions for my five phases
    c1v2 = (xnames=="1")-(xnames=="2") ,
    c1v3 = (xnames=="1")-(xnames=="3") ,
    c1v4 = (xnames=="1")-(xnames=="4") ,
    c1v5 = (xnames=="1")-(xnames=="5") ,
    c2v3 = (xnames=="2")-(xnames=="3") ,
    c2v4 = (xnames=="2")-(xnames=="4") ,
    c2v5 = (xnames=="2")-(xnames=="5") ,
    c3v4 = (xnames=="3")-(xnames=="4") ,
    c3v5 = (xnames=="3")-(xnames=="5") ,
    c4v5 = (xnames=="4")-(xnames=="5")
 
  )
}

# Specify the data in a form that is compatible with BRugs model, as a list:
ySDorig = sd(y)
yMorig = mean(y)
z = ( y - yMorig ) / ySDorig
dataList = list(
  y = z ,
  x = x ,
  Ntotal = Ntotal ,
  NxLvl = NxLvl
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

theData = data.frame( y=dataList$y , x=factor(x,labels=xnames) )
a0 = mean( theData$y )
a = aggregate( theData$y , list( theData$x ) , mean )[,2] - a0
ssw = aggregate( theData$y , list( theData$x ) ,
                function(x){var(x)*(length(x)-1)} )[,2]
sp = sqrt( sum( ssw ) / length( theData$y ) )
initsList = list( a0 = a0 , a = a , sigma = sp , aSD = sd(a) )

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c( "a0" ,  "a" , "b0" , "b" , "sigma" , "aSD" )  
adaptSteps = 1000              # Number of steps to "tune" the samplers.
burnInSteps = 5000            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=100000           # Total number of steps in chains to save.
thinSteps=1                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# Create, initialize, and adapt the model:
jagsModel = jags.model( "model.txt" , data=dataList , inits=initsList , 
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                            n.iter=nPerChain , thin=thinSteps )
# resulting codaSamples object has these indices: 
#   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

checkConvergence = FALSE
if ( checkConvergence ) {
  openGraph(width=7,height=7)
  autocorr.plot( codaSamples[[1]][,c("sigma","aSD")] )
  show( gelman.diag( codaSamples ) )
  effectiveChainLength = effectiveSize( codaSamples ) 
  show( effectiveChainLength )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )

# Extract parameter values
sigmaSample = mcmcChain[,"sigma"] * ySDorig
aSDSample = mcmcChain[,"aSD"] * ySDorig
# Extract b values:
b0Sample = mcmcChain[, "b0" ]
chainLength = length(b0Sample)
bSample = array( 0 , dim=c( dataList$NxLvl , chainLength ) )
for ( xidx in 1:dataList$NxLvl ) {
   bSample[xidx,] = mcmcChain[, paste("b[",xidx,"]",sep="") ]
}
# Convert from standardized b values to original scale b values:
b0Sample = b0Sample * ySDorig + yMorig
bSample = bSample * ySDorig

# plot the SDs:
openGraph(width=7,height=7)
layout( matrix(1:2,nrow=2) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
histInfo = plotPost( sigmaSample , xlab="sigma" , main="Cell SD" , showMode=T )
histInfo = plotPost( aSDSample , xlab="aSD" , main="a SD" , showMode=T )
#saveGraph( file=paste(fileNameRoot,"SD",sep="") , type="eps" )
# Plot b values:
openGraph(dataList$NxLvl*2.75,2.5)
layout( matrix( 1:dataList$NxLvl , nrow=1 ) )
par( mar=c(3,1,2.5,0) , mgp=c(2,0.7,0) )
for ( xidx in 1:dataList$NxLvl ) {
    histInfo = plotPost( bSample[xidx,] ,
              xlab=bquote(beta*1[.(xidx)]) ,
              main=paste("x:",xnames[xidx])  )
}
#saveGraph( file=paste(fileNameRoot,"b",sep="") , type="eps" )

# Display contrast analyses
nContrasts = length( contrastList )
if ( nContrasts > 0 ) {
   nPlotPerRow = 5
   nPlotRow = ceiling(nContrasts/nPlotPerRow)
   nPlotCol = ceiling(nContrasts/nPlotRow)
   openGraph(3.75*nPlotCol,2.5*nPlotRow)
   layout( matrix(1:(nPlotRow*nPlotCol),nrow=nPlotRow,ncol=nPlotCol,byrow=T) )
   par( mar=c(4,0.5,2.5,0.5) , mgp=c(2,0.7,0) )
   for ( cIdx in 1:nContrasts ) {
       contrast = matrix( contrastList[[cIdx]],nrow=1) # make it a row matrix
       incIdx = contrast!=0
       histInfo = plotPost( contrast %*% bSample , compVal=0 , 
                xlab=paste( round(contrast[incIdx],2) , xnames[incIdx] ,
                            c(rep("+",sum(incIdx)-1),"") , collapse=" " ) ,
                cex.lab = 1.0 ,
                main=paste( "X Contrast:", names(contrastList)[cIdx] ) )
   }
   #saveGraph( file=paste(fileNameRoot,"xContrasts",sep="") , type="eps" )
}

# Display data with posterior predictive distributions
openGraph(width=1.5*NxLvl,height=5)
plot(0,0, 
     xlim=c(0.2,NxLvl+0.1) , xlab="X" , 
     xaxt="n" ,
     ylim=c(min(y)-0.2*(max(y)-min(y)),max(y)+0.2*(max(y)-min(y))) , ylab="Y" ,
     main="Data with Posterior Predictive Distrib.")
axis( 1 , at=1:NxLvl , lab=xnames )
for ( j in 1:NxLvl ) {
  yVals = y[x==j]
  points( rep(j,length(yVals))+runif(length(yVals),-0.03,0.03) , 
          yVals , pch=20 , cex=1.5 , col="red" )
  chainSub = round(seq(1,chainLength,length=20))
  for ( chnIdx in chainSub ) {
    m = b0Sample[chnIdx] + bSample[j,chnIdx]
    s = sigmaSample[chnIdx]
    yl = m-1.96*s
    yh = m+1.96*s
    ycomb=seq(yl,yh,length=201)
    ynorm = dnorm(ycomb,mean=m,sd=s)
    ynorm = 0.75*ynorm/max(ynorm)
    lines( j-ynorm , ycomb , col="skyblue" ) # col=chnIdx )
  }
}
#saveGraph(file=paste(fileNameRoot,"PostPred",sep=""),type="eps")

#==============================================================================
# Do NHST ANOVA and t tests:

theData = data.frame( y=y , x=factor(x,labels=xnames) )
aovresult = aov( y ~ x , data = theData ) # NHST ANOVA
cat("\n------------------------------------------------------------------\n\n")
print( summary( aovresult ) )
cat("\n------------------------------------------------------------------\n\n")
print( model.tables( aovresult , "means" ) , digits=4 )
openGraph(width=7,height=7)
boxplot( y ~ x , data = theData , cex.axis=1.25 , ylab="Y")
cat("\n------------------------------------------------------------------\n\n")
print( TukeyHSD( aovresult , "x" , ordered = FALSE ) )
openGraph(width=7,height=7)
plot( TukeyHSD( aovresult , "x" ) )
if ( T ) {
  for ( xIdx1 in 1:(NxLvl-1) ) {
    for ( xIdx2 in (xIdx1+1):NxLvl ) {
      cat("\n----------------------------------------------------------\n\n")
      cat( "xIdx1 = " , xIdx1 , ", xIdx2 = " , xIdx2 ,
           ", M2-M1 = " , mean(y[x==xIdx2])-mean(y[x==xIdx1]) , "\n" )
      print( t.test( y[x==xIdx2] , y[x==xIdx1] , var.equal=T ) ) # t test
    }
  }
}
cat("\n------------------------------------------------------------------\n\n")

#==============================================================================
