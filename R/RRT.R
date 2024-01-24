

#   based on:
#   aces-dev-master\transforms\ctl\rrt\RRT.ctl



#   "Glow" module constants
RRT_GLOW_GAIN   = 0.05
RRT_GLOW_MID    = 0.08

#   Red modifier constants
RRT_RED_SCALE   = 0.82
RRT_RED_PIVOT   = 0.03
RRT_RED_HUE     = 0
RRT_RED_WIDTH   = 135

#   Desaturation contants
RRT_SAT_FACTOR  = 0.96
RRT_SAT_MAT     = NULL  #;calc_sat_adjust_matrix( RRT_SAT_FACTOR, p.AP1_RGB2Y )
RRT_SAT_MAT_inv = NULL  #;base::solve( RRT_SAT_MAT )

ODT_SAT_FACTOR  = 0.93
ODT_SAT_MAT     = NULL
ODT_SAT_MAT_inv = NULL


RRT.TF          = NULL  # definitely exported

RRTsweetener.TF = NULL  # ?
RedModifier.TF  = NULL  # ?
Glow.TF         = NULL  # ?

#   DCDM.OETF   = NULL
DCDM.EOTF   = NULL
    


    
#   this should be called from .onLoad()    
makeRRTplus <- function()
    {
    #domain  = matrix( c(0,2), 2, 3, dimnames=list(NULL,c('ACES.R','ACES.G','ACES.B')) )
    #range   = matrix( c(0,120), 2, 3, dimnames=list(NULL,c('OCES.R','OCES.G','OCES.B')) )    

    RRT.TF  <<- general.RRT()               # TransferFunction( RRT, RRTinv, domain, range, id='RRT.TF' )
    names(RRT.TF$element)   <<- "RRT.TF"    # change name under-the-hood. Should have a method for this one.
        
        
    RRT_SAT_MAT     <<- calc_sat_adjust_matrix( RRT_SAT_FACTOR, p.AP1_RGB2Y )
    RRT_SAT_MAT_inv <<- base::solve( RRT_SAT_MAT )
    
    ODT_SAT_MAT     <<- calc_sat_adjust_matrix( ODT_SAT_FACTOR, p.AP1_RGB2Y )
    ODT_SAT_MAT_inv <<- base::solve( ODT_SAT_MAT )
    
    
    
    
    
    #   these are for mostly for developer use
    domain.sweet    = matrix( c(0,2), 2, 3, dimnames=list(NULL,c('ACES.R','ACES.G','ACES.B')) )
    range.sweet     = matrix( c(0,1), 2, 3, dimnames=list(NULL,c('ACES.R','ACES.G','ACES.B')) )    

    RRTsweetener.TF <<- TransferFunction( rrt_sweeteners, inv_rrt_sweeteners, domain.sweet, range.sweet, id='RRTsweetener.TF' )

    RedModifier.TF  <<- TransferFunction( redmod, redmodinv_precise, domain.sweet, range.sweet, id='RedModifier.TF' )

    Glow.TF         <<- TransferFunction( glow, glowinv, domain.sweet, range.sweet, id='Glow.TF' )

    
    #domain  = matrix( c(0,1), 2, 1, dimnames=list(NULL,'OCES') )
    #range   = matrix( c(0,1), 2, 1, dimnames=list(NULL,'display') )
    #DCDM.OETF       <<- TransferFunction( dcdm_encode, dcdm_decode, domain, range, id="DCDM.OETF" )
    
    domain  = matrix( c(0,1), 2, 1, dimnames=list(NULL,'signal') )
    range   = matrix( c(0,1), 2, 1, dimnames=list(NULL,'display linear') )
    DCDM.EOTF       <<- TransferFunction( dcdm_decode, dcdm_encode, domain, range, id="DCDM.EOTF" )

    return(TRUE)
    }


#   parameterized RRT

general.RRT  <-  function( glowmod="1.1", redmod="1.1" )      # "1.1+pinv"
    {
    if( is.null(glowmod)  ||  is.na(glowmod) )
        glowmodifier = FALSE
    else if( is.logical(glowmod) )
        glowmodifier = glowmod
    else if( is.character(glowmod) )
        glowmodifier = TRUE         #   only one version is known
    else
        {
        log_string( ERROR, "glowmod='%s' is invalid.", as.character(glowmod)[1] )
        return(NULL)
        }

        
    if( is.null(redmod)  ||  is.na(redmod) )
        redmodifier = FALSE
    else if( is.logical(redmod) )
        {
        redmodifier = redmod
        preciseinv  = FALSE
        }
    else if( is.character(redmod) )
        {
        redmodifier = TRUE         #   only one version is known
        preciseinv  = grepl( "pinv$", redmod )
        }
    else
        {
        log_string( ERROR, "redmod='%s' is invalid.", as.character(redmod)[1] )
        return(NULL)
        }

    #   Reference Rendering Transform (RRT)
    #
    #   Input  is ACES      AP0
    #   Output is OCES      AP0 as well !

    RRT <- function( aces )
        {
        #--- Glow module     AP0 ---#
        if( glowmodifier )  aces = glow( aces )
        
        #--- Red modifier    AP0 ---#
        if( redmodifier )   aces = redmod( aces )

        #---   AP0 to AP1  (AP1 is RGB rendering space) both are scene-linear ---# 
        # aces    = pmax( aces, 0 )      #   avoids saturated negative colors from becoming positive in the matrix
        rgbPost = as.numeric( p.AP0_2_AP1_MAT %*% aces  )   #; cat( 'rgbPost=', rgbPost, '\n' )

        #--- Global desaturation,   AP1 --- //
        rgbPost  = as.numeric( RRT_SAT_MAT %*% rgbPost )    #; cat( 'rgbPost=', rgbPost, '\n' )

        #--- Apply the tonescale independently in rendering-space RGB     AP1 ---#
        rgbPost = base::sapply( rgbPost, segmented_spline_fwd, C=p.RRT_PARAMS )     #; cat( 'rgbPost=', rgbPost, '\n' )
        #rgbPost[1] = segmented_spline_fwd( rgbPost[1] )
        #rgbPost[2] = segmented_spline_fwd( rgbPost[2] )
        #rgbPost[3] = segmented_spline_fwd( rgbPost[3] )

        #--- RGB rendering space to OCES,  AP0 ---#
        oces    = tcrossprod( rgbPost, p.AP1_2_AP0_MAT )   #; cat( 'oces=', oces, '\n' )

        return( oces )  # still in AP0
        }

    #   Reference Rendering Transform (RRT) inverse
    #
    #   input is OCES       AP0
    #   output is ACES      AP0 as well !

    RRTinv <- function( oces )
        {
        #--- OCES to AP1  ---#
        rgbPost = as.numeric( tcrossprod( oces, p.AP0_2_AP1_MAT ) )

        #--- Apply the tonescale independently in rendering-space RGB     AP1 ---#
        rgbPost = base::sapply( rgbPost, segmented_spline_rev, C=p.RRT_PARAMS )
        #rgbPost[1] = segmented_spline_fwd( rgbPost[1] )
        #rgbPost[2] = segmented_spline_fwd( rgbPost[2] )
        #rgbPost[3] = segmented_spline_fwd( rgbPost[3] )
        
        #--- Global desaturation,   AP1 --- //
        rgbPost = as.numeric( RRT_SAT_MAT_inv %*% rgbPost )

        #---   AP1 to AP0  (AP1 is RGB rendering space) both are scene-linear ---# 
        # aces    = pmax( aces, 0 )      #   avoids saturated negative colors from becoming positive in the matrix
        aces    = as.numeric( p.AP1_2_AP0_MAT %*% rgbPost  )
        
        
        if( redmodifier )
            {
            #--- Red modifier inverse   AP0 ---#
            if( preciseinv )
                aces = redmodinv_precise( aces )
            else
                aces = redmodinv( aces )    # rough approximation
            }
            
            
        #--- Glow module inverse    AP0 ---#
        if( glowmodifier )  aces = glowinv( aces )

        return( aces )  # still in AP0
        }
     
    domain  = matrix( c(0,47000), 2, 3, dimnames=list(NULL,c('ACES.R','ACES.G','ACES.B')) )
    range   = matrix( c(0,10000), 2, 3, dimnames=list(NULL,c('OCES.R','OCES.G','OCES.B')) )    

    TransferFunction( fun=RRT, funinv=RRTinv, domain, range, id=sigfunction() )
    }

#   ------- Glow module functions  --------------  #
glow_fwd    <- function( ycIn, glowGainIn, glowMid )
    {
    if (ycIn <= (2/3) * glowMid )
        glowGainOut = glowGainIn
    else if ( ycIn >= 2 * glowMid )
        glowGainOut = 0
    else
        glowGainOut = glowGainIn * (glowMid / ycIn - 1/2)

    return( glowGainOut )
    }

glow_inv    <- function( ycOut, glowGainIn, glowMid )
    {
    if( ycOut <= ((1 + glowGainIn) * (2/3) * glowMid) ) 
        glowGainOut = -glowGainIn / (1 + glowGainIn)
    else if( ycOut >= (2 * glowMid) )
        glowGainOut = 0
    else 
        glowGainOut = glowGainIn * (glowMid / ycOut - 1/2) / (glowGainIn / 2 - 1)
 
    return( glowGainOut )
    }

    
sigmoid_shaper <- function(x)
    {
    #   Sigmoid function in the range 0 to 1 spanning -2 to +2.
    t = max(1 - abs(x/2), 0)
    y = 1 + sign(x) * (1 - t*t)

    return( y / 2 )
    }
    
    

#-------     Red modifier functions    --------------  #
#   x   input scalar
#   w   full base width of the shaper function (in degrees)
cubic_basis_shaper <- function( x, w )
    {
    M = c( c( -1,  3, -3,  1 ),
           c(  3, -6,  3,  0 ),
           c( -3,  0,  3,  0 ),
           c(  1,  4,  1,  0 ) )
           
    M   = matrix( M/6, 4, 4, byrow=TRUE )
  
    knots   = c( -w/2, -w/4, 0, w/4, w/2 )
  
    if( x <= knots[1]  ||  knots[5] <= x )    return(0)
  

    knot_coord = (x - knots[1]) * 4/w 
    j   = floor(knot_coord)      
    
    # j must be 0,1,2, or 3; but consider float arithmetic (is w/2+w/2==w ?) and test for 4 anyway
    if( j == 4 )    return(0)
    
    t   = knot_coord - j
      
    monomials   = c( t*t*t, t*t, t, 1)
    
    y   = sum( monomials * M[ ,4-j] )
  
    return( y * 3/2 )
    }
    
    
#   translate (hue - centerH) to interval [-180,180)
center_hue <- function( hue, centerH )
    {
    hueCentered = hue - centerH
    
    return( ((hueCentered+180) %% 360) - 180 )
    }

#   translate (hueCentered + centerH) to interval [0,360)    
uncenter_hue <- function( hueCentered, centerH )
    {
    return( (hueCentered + centerH) %% 360 )
    }

#   aces    AP0    
glow    <- function( aces )
    {
    #   --- Glow module ---  This only seems to affect colors near black 
    saturation  = rgb_2_saturation( aces )
    ycIn        = rgb_2_yc( aces )
    s           = sigmoid_shaper( (saturation - 0.4) / 0.2)
    addedGlow   = 1 + glow_fwd( ycIn, RRT_GLOW_GAIN * s, RRT_GLOW_MID )

    #if( addedGlow != 1 )    
    #    {
    #    cat( aces, ' ', addedGlow, '\n' )
    #    }
        
        
    return( addedGlow * aces )
    }
    
#   aces    AP0 RGB       
glowinv    <- function( aces )
    {    
    #   --- Glow module ---     
    saturation  = rgb_2_saturation( aces )
    ycOut       = rgb_2_yc( aces)
    s           = sigmoid_shaper( (saturation - 0.4) / 0.2 )
    reducedGlow = 1 + glow_inv( ycOut, RRT_GLOW_GAIN * s, RRT_GLOW_MID )

    return( reducedGlow * aces )
    }
    

redmod <- function( aces )
    {
    #--- Red modifier ---#
    hue         = rgb_2_hue( aces )                                 #; print(hue)
    centeredHue = center_hue( hue, RRT_RED_HUE )                    #; print(centeredHue)
    hueWeight   = cubic_basis_shaper( centeredHue, RRT_RED_WIDTH )  #; print(hueWeight)
    
    saturation  = rgb_2_saturation( aces )                          #; print(saturation)
    
    aces[1] = aces[1] + hueWeight * saturation * (RRT_RED_PIVOT - aces[1]) * (1 - RRT_RED_SCALE)
    
    return(aces)
    }
    
    
redmodinv <- function( aces )
    {    
    #--- Red modifier inverse ---#
    hue         = rgb_2_hue( aces )
    centeredHue = center_hue( hue, RRT_RED_HUE ) 
    hueWeight   = cubic_basis_shaper( centeredHue, RRT_RED_WIDTH ) 

    if (centeredHue < 0) { #min_f3(aces) = aces[1] (i.e. magenta-red)
      minChan = aces[2] # green
    } else { # min_f3(aces) = aces[2] (i.e. yellow-red)
      minChan = aces[3] # blue
    }

    a   = hueWeight * (1 - RRT_RED_SCALE) - 1
    b   = aces[1] - hueWeight * (RRT_RED_PIVOT + minChan) * (1 - RRT_RED_SCALE)
    c   = hueWeight * RRT_RED_PIVOT * minChan * (1 - RRT_RED_SCALE)

    d   =  b * b - 4 * a * c
    if( d < 0 )
        {
        print( d )
        print( aces )
        }
        
    aces[1] = ( -b - sqrt(d) ) / (2*a)

    return(aces)
    }
    
redmodinv_precise <- function( aces )
    {    
    #--- Red modifier inverse, with precise rootfinder  -  stats::uniroot() ---#
    
    if( aces[1] == RRT_RED_PIVOT )  return(aces)    # no change
    
    myfun   <- function( red, rgb )
        {
        acestest    = c( red, rgb[2:3] )
        return( redmod( acestest )[1] - aces[1] )
        }
    
    if( aces[1] < RRT_RED_PIVOT )
        interval    = c(0,RRT_RED_PIVOT)
    else
        interval    = c(0,1.5*aces[1])
    
    #   check endpoint values
    f1  = myfun( interval[1], aces ) 
    if( f1 == 0 )
        {
        aces[1] = interval[1]
        return( aces )
        }
        
    f2  = myfun( interval[2], aces )    #    ; print(f2)
    
    if( 0 < f1 * f2 )    
        {
        #   same sign
        log_string( WARN, "root-finding interval [%g,%g] invalid.  redmod cannot be inverted.", interval[1], interval[2] )
        return( NA_real_ )
        }
    
    res = try( stats::uniroot( myfun, interval=interval, rgb=aces, tol=.Machine$double.eps^0.5 ),  silent=FALSE )
    
    if( inherits(res,"try-error" ) )        #class(res) == "try-error" 
        {
        cat( 'stats::uniroot()  res = ', utils::str(res), '\n', file=stderr() )
        return( NA_real_ )
        }
        
    #   cat( res$iter, '\n' )
        
    aces = c( res$root, aces[2:3] )
            
    return( aces )
    }
    
    
    
#   aces    AP0 RGB           
rrt_sweeteners <- function( aces )
    {    
    #   --- Glow module ---     
    aces = glow( aces )

    #--- Red modifier ---#    
    aces = redmod( aces )
    
    #---   AP0 to AP1  (AP1 is RGB rendering space) ---#
    #   moved these to general.OOTF()
    #aces    = pmax( aces, 0 )   # clamp_f3( aces, 0., HALF_POS_INF);
    #rgbPre  = p.AP0_2_AP1_MAT %*% aces
    #rgbPre  = pmax( rgbPre, 0 ) # clamp_f3( rgbPre, 0., HALF_MAX);
    
    #--- Global desaturation ---#
    #   moved this to general.OOTF()    
    #rgbPre  = as.numeric( RRT_SAT_MAT %*% rgbPre )
    
    return( aces )
    }
    
    
#   aces    AP0 RGB       
inv_rrt_sweeteners <- function( aces )  # aces was named rgbPost
    {
    #--- Global desaturation ---#
    #   moved this to general.OOTF()        
    #rgbPost = RRT_SAT_MAT_inv %*% rgbPost

    #--- AP1 (ACES RGB rendering space) to AP0 ---#    
    #   moved these to general.OOTF()            
    #rgbPost = pmax( rgbPost, 0 )    #clamp_f3( rgbPost, 0., HALF_MAX);
    #aces    = p.AP1_2_AP0_MAT %*% rgbPost
    #aces    = pmax( aces, 0 )       #clamp_f3( aces, 0., HALF_MAX);

    #--- Red modifier inverse ---#
    aces = redmodinv_precise( aces )

    #--- Glow module inverse ---#
    aces = glowinv( aces )
    
    return( aces )
    }
    
    