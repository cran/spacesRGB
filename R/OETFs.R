
#   the *special* OETFs are:
#       sRGB
#       BT.709
#       BT.2020
#       ProPhotoRGB
#       240M                (ANSI/SMPTE 240M)

    
#############     sRGB  OETF and OETFinv [0,1]  <-->  [0,1]   ###############    
#   maps [0,1] to [0,1].  
#   it is OK if input is a matrix, and then the return value is a matrix of the same shape
OETF_sRGB <- function( lin )
    {
    out = ifelse( lin <= 0.0031308,    12.92 * lin,  (1055 * lin^(1/2.4) - 55 ) / 1000 )
    return( out )
    }

OETFinv_sRGB <- function( sig )
    {
    out = ifelse( sig <= 12.92*0.0031308,   sig / 12.92,   ( (1000*sig + 55)/(1055) ) ^ 2.4  )
    return( out )
    }    
    

#############     BT.709  OETF and OETFinv [0,1]  <-->  [0,1]   ###############    
#   maps [0,1] to [0,1].  
#   it is OK if input is a matrix, and then the return value is a matrix of the same shape
#   NOTE: 4.5 has been changed to 4.513775 so that both functions are monotone
OETF_BT.709 <- function( lin )
    {
    ifelse( lin < 0.018,  4.513775*lin,  1.099*(lin^0.45 - 1) + 1 )
    }

OETFinv_BT.709 <- function( sig )
    {
    ifelse( sig < 4.513775*0.018,  sig/4.513775,  ( (sig+0.099)/1.099 )^(1/0.45)  )
    }    
    
    
    
#############     BT.2020  OETF and OETFinv [0,1]  <-->  [0,1]   ###############    
#   maps [0,1] to [0,1].  
#   it is OK if input is a matrix, and then the return value is a matrix of the same shape
#   NOTE: 4.5 has been changed to 4.513775 so that both functions are monotone
OETF_BT.2020 <- function( lin )
    {
    alpha   = 1.09929682680944
    beta    = 0.018053968510807
    ifelse( lin < beta,  4.5*lin,  alpha*(lin^0.45 - 1) + 1 )
    }

OETFinv_BT.2020 <- function( sig )
    {
    alpha   = 1.09929682680944
    beta    = 0.018053968510807    
    ifelse( sig < 4.5*beta,  sig/4.5,  ( ((sig-1) + alpha)/alpha )^(1/0.45)  )
    }    



#############     ProPhotoRGB  OETF and OETFinv [0,1]  <-->  [0,1]   ###############    
#   maps [0,1] to [0,1].  
#   it is OK if input is a matrix, and then the return value is a matrix of the same shape
OETF_ProPhotoRGB <- function( lin )
    {
    out = ifelse( lin <= 1/512,  16 * lin,  lin^(1/1.8) )
    return( out )
    }

OETFinv_ProPhotoRGB <- function( sig )
    {
    out = ifelse( sig <= 1/32,   sig/16,  sig ^ 1.8  )
    return( out )
    }    
    
    
#############     240M  OETF and OETFinv [0,1]  <-->  [0,1]   ###############    
#   maps [0,1] to [0,1].  
#   the slope was changed from 4 to 4.002588 in order to pass the round-trip tests
#   it is OK if input is a matrix, and then the return value is a matrix of the same shape
OETF_240M <- function( lin )
    {
    ifelse( lin < 0.0228,  4.002588*lin, 1.1115*(lin^0.45 - 1) + 1 )     # (11115*lin^0.45 - 1115)/10000 )
    }

OETFinv_240M <- function( sig )
    {
    ifelse( sig < 4.002588 * 0.0228,   sig/4.002588,  ( (sig+0.1115)/1.1115 )^(1/0.45)  )
    }    
        
        
#   TF      a transfer function to be validated
#   things checked are
#       *) maps 0 to 0 and 1 to 1
#       *) is strictly monotone        
#       *) preserves dimensions
validTF  <-  function( TF ) 
    {
    #   check that 0->0 and 1->1
    endpoint    = c(0,1)
    ok  = all( TF(endpoint) == endpoint )
    if( ! ok )
        {
        log.string( ERROR, "Transfer Function does not map 0->0 and 1->1." )
        return(FALSE)
        }    
    
    x   = seq( 0, 1, by=1/4096 )
    y   = TF( x )
    ok  = all( 0 < diff(y) )
    if( ! ok )
        {
        log.string( ERROR, "Transfer Function is not strictly monotone." )
        return(FALSE)
        }        
        

    dim(x)   = c( 17, 241 )
    y   = TF( x )    
    ok  = all( dim(y) == dim(x) )
    if( ! ok )
        {
        log.string( ERROR, "Transfer Function does not preserve dimensions." )
        return(FALSE)
        }        

    return(TRUE)
    }
    
    
#   check that TF1 and TF2 are inverses, on the interval [0,1]
validTF_pair  <-  function( TF1, TF2, digits=5 )        
    {    
    tol = 0.5 * 10^(-digits)
    
    x       = seq( 0, 1, by=1/4096 )
    
    x.back  = TF1( TF2(x) )
    delta   = max( abs( x - x.back ) )  #; print(delta)
    if( tol <= delta )
        {
        log.string( ERROR, "Transfer Function pair are not inverses, to %d digits.  TF1(TF2(x))", digits )
        return(FALSE)
        }       
    
    x.back  = TF2( TF1(x) )
    delta   = max( abs( x - x.back ) )
    if( tol <= delta )
        {
        log.string( ERROR, "Transfer Function pair are not inverses, to %d digits.  TF2(TF1(x))", digits )
        return(FALSE)
        }       

    return( TRUE )
    }
    
    
#   TF      a transfer function that has already been validated
makeInverseTF  <-  function( TF )
    {
    #   estimate the rough gamma of TF
    y       = TF(0.5)
    gamma   = log(y) / log(0.5) #; print(gamma)
    
    x   = seq( 0, 1, by=1/512 ) ^ (1/gamma)
    y   = TF(x)
    
    out = splinefun( y, x, method='monoH.FC' )
    
    return( out )
    }
    
#   TF      a transfer function that has already been validated
#   returns optimal gamma in the L1 norm
fittedGammaL1 <- function( TF )
    {
    x   = seq( 0, 1, by=1/512 )
    ytf = TF( x )
    
    myfun <- function( gamma )
        {
        return( sum( abs(x^gamma - ytf) ) )
        }    
    
    #   make initial estimate of gamma
    gamma   = log( TF(0.5) ) / log(0.5)    
        
    res = optimize( myfun, lower=0.5*gamma, upper=2*gamma ) #; print( str(res) )
        
    #   log.string( INFO, "LM polished gamma=%g to %g in %d iterations.", gamma, res$par, res$niter )

    return( res$minimum )
    }
    
    
    
#   TF      a transfer function that has already been validated
#   returns optimal gamma in the least-squares norm
#fittedGammaL2 <- function( TF )
#    {
#    if( ! requireNamespace( 'minpack.lm', quietly=TRUE ) )
#        {
#        log.string( ERROR, "Best-fit gamma canot be computed, because package 'minpack.lm' cannot be loaded." )
#        return(  NA_real_ )
#        }
#    
#    x   = seq( 0, 1, by=1/512 )
#    ytf = TF( x )
#    
#    myresid <- function( gamma )
#        {
#        return( x^gamma - ytf )
#        }    
#    
#    #   make initial estimate of gamma
#    gamma   = log( TF(0.5) ) / log(0.5)    
#        
#    #  comment this out to avoid a WARNING
#    res = minpack.lm::nls.lm( gamma, lower=0.1, upper=10, fn=myresid ) #; print( str(res) )
#    
#    ok  = (1 <= res$info) && (res$info <= 4)
#    if( ! ok )
#        {
#        log.string( WARN, "Levenberg-Marquardt did not converge. info=%d",  res$info )
#        return( NA_real_ )
#        }
#        
    #   log.string( INFO, "LM polished gamma=%g to %g in %d iterations.", gamma, res$par, res$niter )
#
#    return( res$par )
#    }
    
    
    