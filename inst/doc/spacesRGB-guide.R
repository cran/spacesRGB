## ----setup, include=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options( width=100 )
library(spacesRGB)

## ----echo=TRUE, message=FALSE, fig.cap='Figure 1.2  The 3 Transfer Functions of sRGB', fig.align="center", fig.width=4, fig.height=4, dev='png'----
theSpace = getRGB('sRGB')
par( omi=c(0,0,0,0), mai=c(0.6,0.55,0.1,0.1) )
plot( theSpace$OETF, main='', ylab='', color='red' )
plot( theSpace$EOTF, add=TRUE, color='blue' ) ; plot( theSpace$OOTF, add=TRUE, color='black' )
legend( 'topleft', legend=c('OETF','EOTF','OOTF'), bty='n', lwd=2, col=c('red','blue','black') )

## ----echo=TRUE, message=TRUE----------------------------------------------------------------------
prim  = matrix( c(0.64,0.33,  0.30,0.60,  0.15,0.06,  0.3127,0.3290), 4, 2, byrow=TRUE )
installRGB( 'HD+2.4', scene=prim, OETF=BT.709.EOTF^-1, EOTF=BT.1886.EOTF(), overwrite=TRUE )  

## ----echo=TRUE, message=FALSE, fig.cap='Figure 2.1  The 3 Transfer Functions of HD+2.4', fig.align="center", fig.width=4, fig.height=4, dev='png'----
theSpace = getRGB('HD+2.4')
par( omi=c(0,0,0,0), mai=c(0.6,0.55,0.1,0.1) )
plot( theSpace$OETF, main='', ylab='', color='red' )
plot( theSpace$EOTF, add=TRUE, color='blue' ) ; plot( theSpace$OOTF, add=TRUE, color='black' )
legend( 'topleft', legend=c('OETF','EOTF','OOTF'), bty='n', lwd=2, col=c('red','blue','black') )

## ----echo=TRUE, message=FALSE, fig.cap='Figure 3.1  Building TransferFunctions from Simpler Parts', fig.align="center", fig.width=4, fig.height=4, dev='png'----
# create the squaring function on [0,1]
squaring  = TransferFunction( function(x) {x*x}, sqrt, domain=c(0,1), range=c(0,1) )
par( omi=c(0,0,0,0), mai=c(0.6,0.55,0.1,0.1) )
plot( squaring, main='', color='red' )
plot( inverse(squaring), add=TRUE, color='blue' )  # inverse is sqrt()
#   in the next line, power.EOTF(2) is also squaring, so comp should be identity
comp = squaring * inverse( power.EOTF(2.0) )
plot( comp, add=TRUE, color='black' )    
transfer( comp, seq(0,1,length.out=11) )  # verify that comp is the identity

## ----echo=TRUE, message=TRUE----------------------------------------------------------------------
OETF = RRT.TF * general.PODT(REC709_PRI,observer=P3D60_PRI['W',],surround='dim') * sRGB.EOTF^-1 
installRGB( 'sRGB_D60sim', scene=AP0_PRI, EOTF=sRGB.EOTF, OETF=OETF )

## ----echo=TRUE, message=TRUE----------------------------------------------------------------------
installRGB( 'Rec.2020_HLG', scene=AP0_PRI, OOTF=general.OOTF(REC2020_PRI,Ymid=15,Ymax=1000),
						EOTF=HLG.OETF()^-1 * HLG.OOTF(Lw=1) )

## ----echo=FALSE, results='asis'-------------------------------------------------------------------
sessionInfo()

