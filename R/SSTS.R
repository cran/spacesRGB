
#
#   Single-Stage Tone-Scale
#

p.HALF_MIN      = 5.96046448e-08
p.HALF_MAX      = 65504

MIN_STOP_SDR    = -6.5      # only used in this file
MAX_STOP_SDR    =  6.5      # only used in this file

MIN_STOP_RRT    = -15       # only used in this file
MAX_STOP_RRT    =  18       # only used in this file

MIN_LUM_SDR     =  0.02     # only used in this file
MAX_LUM_SDR     = 48.0      # only used in this file

MIN_LUM_RRT     = 0.0001    # only used in this file
MAX_LUM_RRT     = 10000     # only used in this file


#   Textbook monomial to basis-function conversion matrix.
# M1  = matrix(  c( 0.5,-1.0,0.5,  -1.0,1.0,0.5,   0.5,0.0,0.0), 3, 3, byrow=TRUE )



init_TsParams   <- function( minLum, maxLum, expShift=0 )
    {
    ok  = minLum < maxLum
    if( ! ok )
        {
        log_level( ERROR, "minLum=%g and maxLum=%g are invalid.", minLum, maxLum )
        return(NULL)
        }
    
    
    MIN_PT  = list( x=lookup_ACESmin(minLum), y=minLum, slope=0 )
    MID_PT  = list( x=0.18, y=4.8, slope=1.55 )
    MAX_PT  = list( x=lookup_ACESmax(maxLum), y=maxLum, slope=0 )
    
    cLow    = init_coefsLow( MIN_PT, MID_PT)
    cHigh   = init_coefsHigh( MID_PT, MAX_PT)
    
    MIN_PT$x = shift( lookup_ACESmin(minLum), expShift )
    MID_PT$x = shift( 0.18, expShift )
    MAX_PT$x = shift( lookup_ACESmax(maxLum), expShift )

    P   = list( Min=MIN_PT, Mid=MID_PT, Max=MAX_PT, coefsLow=c( cLow, cLow[5] ), coefsHigh=c( cHigh, cHigh[5] ) )
         
    return( P )
    }

    
init_TsParams3  <- function( minLum, midLum, maxLum )
    {
    ok  = minLum < midLum  &&  midLum < maxLum
    if( ! ok )
        {
        log_level( ERROR, "minLum=%g and midLum=%g and maxLum=%g are invalid.", minLum, midLum, maxLum )
        return(NULL)
        }
    
    #   NOTE: This is a bit of a hack - probably a more direct way to do this.
    #   TODO: Fix in future version
    PARAMS_DEFAULT = init_TsParams( minLum, maxLum )
    
    expShift = log2( inv_ssts(midLum,PARAMS_DEFAULT) ) - log2(0.18)
    
    PARAMS = init_TsParams( minLum, maxLum , expShift )
    
    return( PARAMS )
    }
    


lookup_ACESmin    <- function( minLum )
    {
    minTable    = matrix( c( log10(MIN_LUM_RRT), MIN_STOP_RRT, log10(MIN_LUM_SDR), MIN_STOP_SDR ), 2, 2, byrow=TRUE )

    return( 0.18 * 2 ^ interpolate1D( minTable, log10(minLum) ) )
    }
    
lookup_ACESmax    <- function( maxLum )
    {
    maxTable    = matrix( c( log10(MAX_LUM_SDR), MAX_STOP_SDR, log10(MAX_LUM_RRT), MAX_STOP_RRT ), 2, 2, byrow=TRUE )

    return( 0.18 * 2 ^ interpolate1D( maxTable, log10(maxLum) ) )
    }    
    
shift   <- function( x, expShift )
    {
    return( x * 2^(-expShift) )     #2 ^ (log2(x) - expShift) )
    }

init_coefsLow   <- function( TsPointLow, TsPointMid )
    {
    coefsLow    = numeric(5)

    knotIncLow = (log10(TsPointMid$x) - log10(TsPointLow$x)) / 3 

    #   Determine two lowest coefficients (straddling minPt)
    coefsLow[1] = (TsPointLow$slope * (log10(TsPointLow$x)-0.5*knotIncLow)) + ( log10(TsPointLow$y) - TsPointLow$slope * log10(TsPointLow$x) )
    coefsLow[2] = (TsPointLow$slope * (log10(TsPointLow$x)+0.5*knotIncLow)) + ( log10(TsPointLow$y) - TsPointLow$slope * log10(TsPointLow$x) )
    
    #   NOTE: if slope=0, then the above becomes just 
        #   coefsLow[0] = log10(TsPointLow.y);
        #   coefsLow[1] = log10(TsPointLow.y);
    #   leaving it as a variable for now in case we decide we need non-zero slope extensions

    #   Determine two highest coefficients (straddling midPt)
    coefsLow[4] = (TsPointMid$slope * (log10(TsPointMid$x)-0.5*knotIncLow)) + ( log10(TsPointMid$y) - TsPointMid$slope * log10(TsPointMid$x) )
    coefsLow[5] = (TsPointMid$slope * (log10(TsPointMid$x)+0.5*knotIncLow)) + ( log10(TsPointMid$y) - TsPointMid$slope * log10(TsPointMid$x) )
    
    #   Middle coefficient (which defines the "sharpness of the bend") is linearly interpolated
    bendsLow    = matrix( c( MIN_STOP_RRT, 0.18,  MIN_STOP_SDR, 0.35), 2, 2, byrow=TRUE )
    
    pctLow = interpolate1D( bendsLow, log2(TsPointLow$x/0.18) )
    
    #   coefsLow[3] = log10(TsPointLow$y) + pctLow*( log10(TsPointMid$y) - log10(TsPointLow$y) )

    coefsLow[3] = (1 - pctLow)*log10(TsPointLow$y) + pctLow*log10(TsPointMid$y)
    
    return( coefsLow )
    }
    
    
init_coefsHigh  <-  function( TsPointMid, TsPointMax )
    {
    coefsHigh   = numeric(5)

    knotIncHigh = (log10(TsPointMax$x) - log10(TsPointMid$x)) / 3
    
    #   float halfKnotInc = (log10(TsPointMax.x) - log10(TsPointMid.x)) / 6.;

    #   Determine two lowest coefficients (straddling midPt)
    coefsHigh[1] = (TsPointMid$slope * (log10(TsPointMid$x)-0.5*knotIncHigh)) + ( log10(TsPointMid$y) - TsPointMid$slope * log10(TsPointMid$x) )
    coefsHigh[2] = (TsPointMid$slope * (log10(TsPointMid$x)+0.5*knotIncHigh)) + ( log10(TsPointMid$y) - TsPointMid$slope * log10(TsPointMid$x) )

    #   Determine two highest coefficients (straddling maxPt)
    coefsHigh[4] = (TsPointMax$slope * (log10(TsPointMax$x)-0.5*knotIncHigh)) + ( log10(TsPointMax$y) - TsPointMax$slope * log10(TsPointMax$x) )
    coefsHigh[5] = (TsPointMax$slope * (log10(TsPointMax$x)+0.5*knotIncHigh)) + ( log10(TsPointMax$y) - TsPointMax$slope * log10(TsPointMax$x) )
    
    #   NOTE: if slope=0, then the above becomes just
        #   coefsHigh[0] = log10(TsPointHigh.y);
        #   coefsHigh[1] = log10(TsPointHigh.y);
    #   leaving it as a variable for now in case we decide we need non-zero slope extensions
    
    #   Middle coefficient (which defines the "sharpness of the bend") is linearly interpolated
    bendsHigh   = matrix( c(MAX_STOP_SDR, 0.89, MAX_STOP_RRT, 0.90 ), 2, 2, byrow=TRUE )
    
    pctHigh = interpolate1D( bendsHigh, log2(TsPointMax$x/0.18) )
    
    #   coefsHigh[2] = log10(TsPointMid.y) + pctHigh*(log10(TsPointMax.y)-log10(TsPointMid.y));
    
    coefsHigh[3] =  (1 - pctHigh)*log10(TsPointMid$y) + pctHigh*log10(TsPointMax$y)
    
    return( coefsHigh )
    }

    
    
#   x   a single number
#   C   TsParams    

ssts    <- function( x, C )
    {
    if( is.na(x) )  return(NA_real_)
    
    N_KNOTS_LOW     = 4
    N_KNOTS_HIGH    = 4

    #   Check for negatives or zero before taking the log. If negative or zero, set to p.HALF_MIN.
    logx = log10( max(x,p.HALF_MIN) ) #; print(logx) ; print( str(C) )

    if ( logx <= log10(C$Min$x) )
        { 
        logy = logx * C$Min$slope + ( log10(C$Min$y) - C$Min$slope * log10(C$Min$x) )
        } 
    else if( log10(C$Min$x) < logx  &&  logx < log10(C$Mid$x)  )
        {
        knot_coord = (N_KNOTS_LOW-1) * (logx-log10(C$Min$x))/(log10(C$Mid$x)-log10(C$Min$x))
        
        j = floor(knot_coord)
        t = knot_coord - j
        j = j+1

        cf =    C$coefsLow[ j:(j+2) ]

        monomials   = c( t*t, t, 1 )
        logy    = sum( monomials * (cf %*% p.M_monomial_to_basis) )
        } 
    else if( ( log10(C$Mid$x) <= logx ) && ( logx < log10(C$Max$x) ) ) 
        {
        knot_coord  = (N_KNOTS_HIGH-1) * (logx-log10(C$Mid$x))/(log10(C$Max$x)-log10(C$Mid$x));
        j   = floor(knot_coord)
        t   = knot_coord - j
        j   = j+1

        cf  =   C$coefsHigh[ j:(j+2) ]

        monomials   = c( t*t, t, 1 )
        logy = sum( monomials * (cf %*% p.M_monomial_to_basis) )
        } 
    else 
        { # if ( logIn >= log10(C.Max.x) ) { 
        logy = logx * C$Max$slope + ( log10(C$Max$y) - C$Max$slope * log10(C$Max$x) ) 
        }

    return( 10^logy )
    }

    
#   y   a single number
#   C   TsParams

inv_ssts    <- function( y, C )
    {  
    if( is.na(y) )  return(NA_real_)
    
    N_KNOTS_LOW     = 4 
    N_KNOTS_HIGH    = 4    

    KNOT_INC_LOW    = (log10(C$Mid$x) - log10(C$Min$x)) / (N_KNOTS_LOW - 1)
    KNOT_INC_HIGH   = (log10(C$Max$x) - log10(C$Mid$x)) / (N_KNOTS_HIGH - 1)

    #   KNOT_Y is luminance of the spline at each knot
    KNOT_Y_LOW  = (C$coefsLow[1:4] + C$coefsLow[2:5]) / 2

    KNOT_Y_HIGH = (C$coefsHigh[1:4] + C$coefsHigh[2:5]) / 2

    logy    = log10( max(y,1e-10) )
    
    #   cf  = numeric(3)
        
    if( logy <= log10(C$Min$y) )
        {
        logx = log10(C$Min$x)
        } 
    else if( (logy > log10(C$Min$y)) && (logy <= log10(C$Mid$y)) )
        {
        if( logy > KNOT_Y_LOW[1] && logy <= KNOT_Y_LOW[2])
            {
            cf  = C$coefsLow[1:3];  j = 0
            } 
        else if ( logy > KNOT_Y_LOW[2] && logy <= KNOT_Y_LOW[3] )
            {
            cf  = C$coefsLow[2:4];  j = 1
            }
        else if ( logy > KNOT_Y_LOW[3] && logy <= KNOT_Y_LOW[4] )
            {
            cf  = C$coefsLow[3:5];  j = 2
            } 

        tmp = cf %*% p.M_monomial_to_basis

        a   = tmp[1]
        b   = tmp[2]
        c   = tmp[3]
        c   = c - logy

        d = sqrt( b*b - 4*a*c )

        t = (2 * c) / (-d -b)

        logx    = log10(C$Min$x) + (j + t) * KNOT_INC_LOW
        } 
    else if( (logy > log10(C$Mid$y)) && (logy < log10(C$Max$y)) )
        {
        if ( logy >= KNOT_Y_HIGH[1] && logy <= KNOT_Y_HIGH[2])
            {
            cf  = C$coefsHigh[1:3] ;  j = 0
            }
        else if( logy > KNOT_Y_HIGH[2] && logy <= KNOT_Y_HIGH[3])
            {
            cf  = C$coefsHigh[2:4] ;  j = 1
            }
        else if( logy > KNOT_Y_HIGH[3] && logy <= KNOT_Y_HIGH[4])
            {
            cf  = C$coefsHigh[3:5] ;  j = 2
            } 

        tmp = cf %*% p.M_monomial_to_basis

        a   = tmp[1]
        b   = tmp[2]
        c   = tmp[3]
        c   = c - logy

        d = sqrt( b*b - 4*a*c )

        t = (2 * c) / ( -d - b)

        logx    = log10(C$Mid$x) + (t + j) * KNOT_INC_HIGH
        }
    else
        { #//if ( logy >= log10(C.Max.y) ) {
        logx    = log10(C$Max$x)
        }

    return( 10 ^ logx )
    }
    

    
#   SSTS.TF is not exported, it is used for testing
SSTS.TF <- function( TsParams )
    {
    ok  = is.list(TsParams)  &&  !is.null(TsParams$Min)  &&  !is.null(TsParams$Mid)  &&  !is.null(TsParams$Max)
    if( ! ok )
        {
        log_level( ERROR, "TsParams='%s' is invalid.", as.character(TsParams)[1] )
        return(NULL)        
        }
    
    domain  = matrix( c(TsParams$Min$x,TsParams$Max$x), 2, 1, dimnames=list(NULL,"AU") )
    
    fun <- function(x)
        {
        out = rep( NA_real_, length(x) )
        
        for( k in 1:length(x) )
            {
            if( is.na( x[k] ) ) next
            
            out[k]  = ssts( x[k], TsParams )
            }
            
        return( out )
        }
        
    funinv <- function(y)
        {
        out = rep( NA_real_, length(y) )
        
        for( k in 1:length(y) )
            {
            if( is.na( y[k] ) ) next
            
            out[k]  = inv_ssts( y[k], TsParams )
            }
            
        return( out )
        }        
    
    yrange  = fun( domain )
    
    range   = matrix( sort(yrange), 2, 1, dimnames=list(NULL,"AU") )
    
    spacesRGB::TransferFunction( fun, funinv, domain, range )
    }
    
    
    
#   float interpolate1D (float table[][2], float p);

interpolate1D <- function( tab, p )
    {
    stats::approx( tab[ ,1], tab[ ,2], p, rule=2 )$y
    }
    