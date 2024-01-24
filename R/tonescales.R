
#   Textbook monomial to basis-function conversion matrix.
M   = matrix(  c( 0.5,-1.0,0.5,  -1.0,1.0,0.5,   0.5,0.0,0.0), 3, 3, byrow=TRUE )


#   this parameter block is suitable for the former segmented_spline_c5_***()
p.RRT_PARAMS = list(
    coefsLow    = c( -4.0000000000, -4.0000000000, -3.1573765773, -0.4852499958, 1.8477324706, 1.8477324706 ),    # coefsLow[6]
    coefsHigh   = c( -0.7185482425, 2.0810307172, 3.6681241237, 4.0000000000, 4.0000000000, 4.0000000000 ),       #     coefsHigh[6]
    minPoint    = list( x=0.18*2^-15,   y=0.0001),  #   minPoint
    midPoint    = list( x=0.18,         y=4.8),     #   midPoint  
    maxPoint    = list( x=0.18*2^18,    y=10000.),  #   maxPoint
    slopeLow    = 0.0,  #   slopeLow
    slopeHigh   = 0.0   #   slopeHigh
    )



segmented_spline_fwd <- function( x, C )
    {
    if( is.na(x) )  return(NA_real_)
    
    N_KNOTS_LOW     = length(C$coefsLow)  - 2  # 4
    N_KNOTS_HIGH    = length(C$coefsHigh) - 2  # 4

    #   Check for negatives or zero before taking the log. If negative or zero, set to HALF_MIN.
    logx = log10( max(x, HALF_MIN ) ) 


    if( logx <= log10(C$minPoint$x) )
        { 
        logy = logx * C$slopeLow + ( log10(C$minPoint$y) - C$slopeLow * log10(C$minPoint$x) )
        } 
    else if( ( logx > log10(C$minPoint$x) ) && ( logx < log10(C$midPoint$x) )) 
        {
        knot_coord = (N_KNOTS_LOW-1) * (logx-log10(C$minPoint$x))/(log10(C$midPoint$x) - log10(C$minPoint$x))
        j = floor(knot_coord)
        t = knot_coord - j

        cf  = C$coefsLow[ (j+1):(j+3) ]

        monomials   = c( t*t, t, 1 )
        logy        = sum( monomials * (cf %*% M) )
        }
    else if (( logx >= log10(C$midPoint$x) ) && ( logx < log10(C$maxPoint$x) )) 
        {
        knot_coord = (N_KNOTS_HIGH-1) * (logx - log10(C$midPoint$x))/(log10(C$maxPoint$x) - log10(C$midPoint$x))
        j = floor( knot_coord )
        t = knot_coord - j 
        cf  = C$coefsHigh[ (j+1):(j+3) ]
        
        monomials   = c( t*t, t, 1 )
        logy        = sum( monomials * (cf %*% M) )
        } 
    else
        { #if ( logIn >= log10(C.maxPoint.x) ) { 
        logy = logx * C$slopeHigh + ( log10(C$maxPoint$y) - C$slopeHigh * log10(C$maxPoint$x) )
        }

    return( 10^logy )
    }


segmented_spline_rev <- function( y, C )
    {    
    if( is.na(y) )  return(NA_real_)
        
    N_KNOTS_LOW     = length(C$coefsLow)  - 2  # 4
    N_KNOTS_HIGH    = length(C$coefsHigh) - 2  # 4

    KNOT_INC_LOW    = (log10(C$midPoint$x) - log10(C$minPoint$x)) / (N_KNOTS_LOW - 1)
    KNOT_INC_HIGH   = (log10(C$maxPoint$x) - log10(C$midPoint$x)) / (N_KNOTS_HIGH - 1)
  
    #   KNOT_Y is luminance of the spline at each knot
    KNOT_Y_LOW  = (C$coefsLow[1:N_KNOTS_LOW] + C$coefsLow[1 + (1:N_KNOTS_LOW)]) / 2

    KNOT_Y_HIGH = (C$coefsHigh[1:N_KNOTS_HIGH] + C$coefsHigh[1 + (1:N_KNOTS_HIGH)]) / 2        

    logy = log10( max(y,1e-10) )


    if( logy <= log10(C$minPoint$y) ) 
        {
        logx = log10(C$minPoint$x)
        } 
    else if( (log10(C$minPoint$y) < logy)  &&  (logy <= log10(C$midPoint$y)) )
        {
        j   = base::findInterval( logy, KNOT_Y_LOW, all.inside=TRUE ) 
        
        #print(KNOT_Y_LOW) ; print(logy) ; print( logy - KNOT_Y_LOW )
        #print(j)
        
        cf  = C$coefsLow[ j:(j+2) ]     #; print(cf)
        j   = j - 1
        
        tmp = cf %*% M

        a = tmp[1]
        b = tmp[2]
        c = tmp[3]
        c = c - logy
        
        d   =  b*b - 4*a*c
        
        if( d < 0 )
            {
            #   this happens, and d is about -1.e-11.  There may be a mismatch between C$minPoint$y  and  C$coefsLow[1]
            #   print(d)
            d = 0
            }

        d = sqrt( d )

        t = (2*c) / (-d - b)

        logx = log10(C$minPoint$x) + (t + j) * KNOT_INC_LOW
        }
    else if( (log10(C$midPoint$y) < logy) && (logy < log10(C$maxPoint$y)) )
        {
        j   = base::findInterval( logy, KNOT_Y_HIGH, all.inside=TRUE )
        
        cf  = C$coefsHigh[ j:(j+2) ]
        j   = j - 1
        
        tmp = cf %*% M

        a = tmp[1]
        b = tmp[2]
        c = tmp[3]
        c = c - logy;

        d   =  b*b - 4*a*c
        
        if( d < 0 )
            {
            #   print(d)
            d = 0
            }
            
        d = sqrt( d )

        t = (2*c) / (-d - b)

        logx = log10(C$midPoint$x) + (t + j) * KNOT_INC_HIGH
        }
    else
        {   #if ( logy >= log10(C.maxPoint.y) ) {
        logx = log10(C$maxPoint$x)
        }
  
    return( 10^logx )
    }

    
#   these parameter blocks are suitable for the former segmented_spline_c9_***()
p.ODT_48nits = list(
    coefsLow    = c( -1.6989700043, -1.6989700043, -1.4779000000, -1.2291000000, -0.8648000000, -0.4480000000, 0.0051800000, 0.4511080334, 0.9113744414, 0.9113744414), # coefsLow[10]
    coefsHigh   = c(0.5154386965, 0.8470437783, 1.1358000000, 1.3802000000, 1.5197000000, 1.5985000000, 1.6467000000, 1.6746091357, 1.6878733390, 1.6878733390 ),  # coefsHigh[10]
    minPoint    = list( x=segmented_spline_fwd( 0.18 * 2^-6.5, C=p.RRT_PARAMS),  y=0.02 ),  #   minPoint
    midPoint    = list( x=segmented_spline_fwd( 0.18 , C=p.RRT_PARAMS),          y=4.8  ),  #   midPoint  
    maxPoint    = list( x=segmented_spline_fwd( 0.18 * 2^6.5, C=p.RRT_PARAMS ),  y=48.0 ),  #   maxPoint
    slopeLow    = 0.0,  #   slopeLow
    slopeHigh   = 0.04  #   slopeHigh
    )

p.ODT_1000nits  = list(
    coefsLow    = c( -4.9706219331, -3.0293780669, -2.1262, -1.5105, -1.0578, -0.4668, 0.11938, 0.7088134201, 1.2911865799, 1.2911865799 ),  #   coefsLow[10]
    coefsHigh   = c(0.8089132070, 1.1910867930, 1.5683, 1.9483, 2.3083, 2.6384, 2.8595, 2.9872608805, 3.0127391195, 3.0127391195 ),   #   coefsHigh[10]
    minPoint    = list( x=segmented_spline_fwd( 0.18 * 2^-12, C=p.RRT_PARAMS), y=0.0001 ),  #   minPoint
    midPoint    = list( x=segmented_spline_fwd( 0.18,   C=p.RRT_PARAMS ),      y=10.0 ),    #   midPoint  
    maxPoint    = list( x=segmented_spline_fwd( 0.18 * 2^10, C=p.RRT_PARAMS),  y=1000.0 ),  #   maxPoint
    slopeLow    = 3.0,  #   slopeLow
    slopeHigh   = 0.06  #   slopeHigh
    )

p.ODT_2000nits  = list(
    coefsLow  = c( -4.9706219331, -3.0293780669, -2.1262, -1.5105, -1.0578, -0.4668, 0.11938, 0.7088134201, 1.2911865799, 1.2911865799 ),  #  coefsLow[10]
    coefsHigh = c( 0.8019952042, 1.1980047958, 1.5943000000, 1.9973000000, 2.3783000000, 2.7684000000, 3.0515000000, 3.2746293562, 3.3274306351, 3.3274306351 ),  # coefsHigh[10]
    minPoint    = list( x=segmented_spline_fwd( 0.18 * 2^-12, C=p.RRT_PARAMS), y=0.0001 ),  #  minPoint
    midPoint    = list( x=segmented_spline_fwd( 0.18,   C=p.RRT_PARAMS ),      y=10.0 ),    #  midPoint  
    maxPoint    = list( x=segmented_spline_fwd( 0.18 * 2^11, C=p.RRT_PARAMS ), y=2000.0 ),  #  maxPoint
    slopeLow    = 3.0,  #   slopeLow
    slopeHigh   = 0.12  #   slopeHigh
    )
    
p.ODT_4000nits  = list(
    coefsLow    = c( -4.9706219331, -3.0293780669, -2.1262, -1.5105, -1.0578, -0.4668, 0.11938, 0.7088134201, 1.2911865799, 1.2911865799 ),   # coefsLow[10]
    coefsHigh   = c( 0.7973186613, 1.2026813387, 1.6093000000, 2.0108000000, 2.4148000000, 2.8179000000, 3.1725000000, 3.5344995451, 3.6696204376, 3.6696204376 ),  # coefsHigh[10]
    minPoint    = list( x=segmented_spline_fwd( 0.18 * 2^-12, C=p.RRT_PARAMS ), y=0.0001 ), # minPoint
    midPoint    = list( x=segmented_spline_fwd( 0.18,   C=p.RRT_PARAMS ),       y=10.0 ),   # midPoint  
    maxPoint    = list( x=segmented_spline_fwd( 0.18 * 2^12, C=p.RRT_PARAMS  ), y=4000.0),  # maxPoint
    slopeLow    = 3.0,  #   slopeLow
    slopeHigh   = 0.3   #   slopeHigh
    )



    


##-----     TransferFunction wrapper   -   useful for testing       --------#

segmented_spline.TF <- function( PARAMS )
    {
    ok  = is.list(PARAMS)  &&  !is.null( PARAMS$minPoint)  &&  !is.null(PARAMS$midPoint)  &&  !is.null(PARAMS$maxPoint)
    if( ! ok )
        {
        log_string( ERROR, "PARAMS='%s' is invalid.", as.character(PARAMS)[1] )
        return(NULL)        
        }
    
    domain  = matrix( c(PARAMS$minPoint$x,PARAMS$maxPoint$x), 2, 1, dimnames=list(NULL,"AU") )
    
    fun <- function(x)
        {
        out = rep( NA_real_, length(x) )
        
        for( k in 1:length(x) )
            {
            out[k]  = segmented_spline_fwd( x[k], C=PARAMS )
            }
            
        return( out )
        }
        
    funinv <- function(y)
        {
        out = rep( NA_real_, length(y) )
        
        for( k in 1:length(y) )
            {
            out[k]  = segmented_spline_rev( y[k], C=PARAMS )
            }
            
        return( out )
        }        
    
    yrange  = fun( domain )
    
    range   = matrix( sort(yrange), 2, 1, dimnames=list(NULL,"AU") )
    
    TransferFunction( fun, funinv, domain, range, id=sigfunction() )
    }
    
    
    