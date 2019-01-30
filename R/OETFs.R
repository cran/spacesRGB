
#   the *special* OETFs are:
#       sRGB
#       BT.709
#       BT.2020
#       ProPhotoRGB
#       240M                (ANSI/SMPTE 240M)


#   The "Gang of 5"
sRGB.EOTF           = NULL
BT.709.EOTF         = NULL
BT.2020.EOTF        = NULL
ProPhotoRGB.EOTF    = NULL
SMPTE.240M.EOTF     = NULL

FullRangeToSMPTE.TF = NULL

#   this should be called from .onLoad()
makeGangOf5   <- function()
    {
    domain              =   matrix( c(0,1), 2, 1, dimnames=list(NULL,"non-linear signal") )
    range               =   matrix( c(0,1), 2, 1, dimnames=list(NULL,"linear display") )
    
    sRGB.EOTF           <<- TransferFunction( OETFinv_sRGB, OETF_sRGB, domain, range, id="sRGB.EOTF" )
        
    BT.709.EOTF         <<- TransferFunction( OETFinv_BT.709, OETF_BT.709, domain, range, id="BT.709.EOTF" )
        
    BT.2020.EOTF        <<- TransferFunction( OETFinv_BT.2020, OETF_BT.2020, domain, range, id="BT.2020.EOTF" )

    ProPhotoRGB.EOTF    <<- TransferFunction( OETFinv_ProPhotoRGB, OETF_ProPhotoRGB, domain, range, id='ProPhotoRGB.EOTF' )

    SMPTE.240M.EOTF     <<- TransferFunction( OETFinv_240M, OETF_240M, domain, range, id='SMPTE.240M.EOTF' )

    #   add a 6th one
    FullRangeToSMPTE.TF   <<- affine.TF( 64/1023, 940/1023 )
    names(FullRangeToSMPTE.TF$element)  <<- "FullRangeToSMPTE.TF"    # change name under-the-hood. Should have a method for this one.    
    }



    
#############     sRGB  OETF and OETFinv [0,1]  <-->  [0,1]   ###############    
#   maps [0,1] to [0,1].  
#   it is OK if input is a matrix, and then the return value is a matrix of the same shape
OETF_sRGB <- function( lin )
    {
    k0  = 0.040448236277108     #   see rootK0() below
    out = ifelse( lin <= k0/12.92,    12.92 * lin,  (1055 * lin^(1/2.4) - 55 ) / 1000 )
    return( out )
    }

OETFinv_sRGB <- function( sig )
    {
    k0  = 0.040448236277108     #   see rootK0() below    
    out = ifelse( sig <= k0,   sig / 12.92,   ( (1000*sig + 55)/(1055) ) ^ 2.4  )
    return( out )
    }    
    
rootK0  <-  function()
    {
    myfun   <- function(x) { 12.92 * ( (1000*x + 55)/1055 )^2.4  -  x }
    
    stats::uniroot( myfun, c(0.04,0.05),    tol = .Machine$double.eps^1, trace=2 )
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
    

        
        
#############   ITU-R BT.1886  EOTF   [0,1] <--> [Lb,Lw]     ###############    

BT.1886.EOTF <- function( gamma=2.4, Lb=0, Lw=1 )
    {
    if( length(gamma)!=1  ||  gamma <= 0 )
        {
        log.string( ERROR, "gamma=%g is invalid.", gamma )
        return(NULL)
        }
        
    ok  = length(Lb)==1  &&  length(Lw)==1  &&  (0 <= Lb)  &&  (Lb < Lw)
    if( ! ok )
        {
        log.string( ERROR, "Lb=%g  and/or  Lw=%g are invalid.", Lb, Lw )
        return(NULL)
        }
            
    Lb1g    = Lb^(1/gamma)
    Lw1g    = Lw^(1/gamma)    
    
    #   these 3 are not necessary
    #denom   = Lw1g - Lb1g
    #a   = denom ^ gamma
    #b   = Lb1g / denom
    
    fun     <- function(V) { ((1-V)*Lb1g + V*Lw1g)^gamma }
    funinv  <- function(L) { (L^(1/gamma) - Lb1g) / (Lw1g - Lb1g) }
    
    domain  =   matrix( c(0,1), 2, 1, dimnames=list(NULL,"display signal") )
    range   =   matrix( c(Lb,Lw), 2, 1, dimnames=list(NULL,"display linear") )
    
    TransferFunction( fun, funinv, domain, range, id=sigfunction() )
    }        
        
        
    
#############     BT.2020  OETF and OETFinv [0,1]  <-->  [0,1]   ###############    
#   maps [0,1] to [0,1].  
#   it is OK if input is a matrix, and then the return value is a matrix of the same shape

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
    

EOTFfromString  <-  function( id )
    {
    ok  = is.character(id)  &&  length(id)==1
    if( ! ok )
        {
        log.string( ERROR, "id='%s' is invalid.", as.character(id)[1] )
        return(NULL)
        }          
    
    id.full = c( 'sRGB', 'BT.709', 'BT.2020', '240M', 'ProPhotoRGB' )

    idx = pmatch( tolower(id), tolower(id.full) )
    if( is.na(idx) )
        {
        log.string( ERROR, "id='%s' is invalid; it may match more than one EOTF.", id )
        return(NULL)
        }          

    listTF  = list( sRGB.EOTF, BT.709.EOTF, BT.2020.EOTF, SMPTE.240M.EOTF, ProPhotoRGB.EOTF )

    return( listTF[[idx]] )
    }
    



    
    
###########     pure gamma classical  [0,1]  <-->  [0,1]    ###################
#   this one is parameterized by gamma, and returns a TransferFunction
#
power.OETF  <-  function( gamma )
    {
    ok  = is.numeric(gamma)  &&  length(gamma)==1  &&  0<gamma
    if( ! ok )
        {
        log.string( ERROR, "gamma = '%s' is invalid.", as.character(gamma) )
        return(NULL)
        }   
        
    domain  = matrix( c(0,1), 2, 1, dimnames=list(NULL,"linear scene") )
    range   = matrix( c(0,1), 2, 1, dimnames=list(NULL,"non-linear signal") )
    
    id  = sprintf( "power.OETF(%g)", gamma )
    
    out = TransferFunction( function(x) {x^(1/gamma)}, function(y) {y^gamma}, domain, range, id=id )
    
    metadata(out)   = list( gamma=gamma )
    
    return( out )
    }
    
power.EOTF  <-  function( gamma )
    {
    ok  = is.numeric(gamma)  &&  length(gamma)==1  &&  0<gamma
    if( ! ok )
        {
        log.string( ERROR, "gamma = '%s' is invalid.", as.character(gamma) )
        return(NULL)
        }   
        
    domain  = matrix( c(0,1), 2, 1, dimnames=list(NULL,"non-linear signal") )
    range   = matrix( c(0,1), 2, 1, dimnames=list(NULL,"linear display") )
    
    id  = sprintf( "power.EOTF(%g)", gamma )
        
    out = TransferFunction( function(x) {x^gamma}, function(y) {y^(1/gamma)}, domain, range, id=id )
    
    metadata(out)   = list( gamma=gamma )
    
    return( out )    
    }
    

power.OOTF  <-  function( gamma )       #   end-to-end OOTF    
    {
    ok  = is.numeric(gamma)  &&  length(gamma)==1  &&  0<gamma
    if( ! ok )
        {
        log.string( ERROR, "gamma = '%s' is invalid.", as.character(gamma) )
        return(NULL)
        }   
        
    domain  = matrix( c(0,1), 2, 1, dimnames=list(NULL,"non-linear signal") )
    range   = matrix( c(0,1), 2, 1, dimnames=list(NULL,"linear display") )
    
    id  = sprintf( "power.OOTF(%g)", gamma )
        
    out = TransferFunction( function(x) {x^gamma}, function(y) {y^(1/gamma)}, domain, range, id=id )
    
    metadata(out)   = list( gamma=gamma )
    
    return( out )    
    }
    
    
    
        
        
###########     affine classical  [0,1]  <-->  [Y0,Y1]    ###################
#   this one is parameterized by Y0 and Y1, and returns a TransferFunction
#
affine.TF  <-  function( y0, y1 )
    {        
    ok  = is.numeric(y0)  &&  is.numeric(y1)  &&  length(y0)==1  &&  length(y1)==1  &&  y0!=y1  # &&  Ymin<Ymax
    if( ! ok )
        {
        log.string( ERROR, "y0='%s' or y1='%s' is invalid, or they are equal.", 
                                as.character(y0)[1], as.character(y1)[1] )
        return(NULL)
        }   
        
    domain  = matrix( c(0,1), 2, 1, dimnames=list(NULL,"AU") )
    range   = matrix( sort(c(y0,y1)), 2, 1, dimnames=list(NULL,"AU") )
    
    fun     <- function(x)  { (1-x)*y0 + x*y1 }
    funinv  <- function(y)  { (y - y0)/(y1 - y0) }
    
    TransferFunction( fun, funinv, domain, range, id=sigfunction() )
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
    
    