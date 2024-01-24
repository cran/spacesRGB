



DIM_SURROUND_GAMMA  = 0.9811

#   returns a TransferFunction

general.OOTF    <- function(
                    display_pri,                #   4x2 matrix
                    Ymin=0.00010,               #   nit
                    Ymid=7.2,                   #   nit
                    Ymax=108,                   #   nit
                    observerWP=NULL,            #   xy chromaticity (2 numbers) of assumed observer whitepoint
                    limiting_pri=NULL,          #   4x2 matrix, NULL means to not limit                    
                    surround='dark',            #   'dark', 'dim', or 'normal'
                    dynrange='SDR',             #   'SDR' or 'HDR'
                    glowmod="1.1",
                    redmod="1.1"
                    )
    {
    #print( sys.status() )
    #theSig = gsub( ' ', '', as.character( as.expression( sys.calls()[[1]] ) ) )
    #theSig = gsub( ' ', '', deparse( sys.calls()[[1]] ) )
    #   print( sigfunction() )
    
    if( is.null(display_pri) )
        {
        log_string( ERROR, "display_pri cannot be NULL." )
        return(NULL)
        }
    
    if( is.null(observerWP) )
        {
        #   get observerWP from display_pri
        observerWP  = display_pri[4, ]
        }

    if( all( observerWP == AP0_PRI[4, ] ) )       # ||  ! is.na(simscale) )
        {
        #   observer whitepoint is not available, or it matches AP0 whitepoint  - so no need to make CAT
        ACEStoObserverCAT = NULL
        ObservertoACESCAT = NULL
        }    
    else
        {
        if( ! requireNamespace( 'spacesXYZ', quietly=TRUE ) )
            {    
            log_string( ERROR, "Cannot adapt from ACES whitepoint to display (assumed observer adapted) whitepoint, because 'spacesXYZ' cannot be loaded." )
            return( NULL )
            }

        #   make Y=1 for both whitepoints
        white.xyY   = rbind( c(AP0_PRI[4, ],1), c(observerWP,1) )
        white.XYZ   = spacesXYZ::XYZfromxyY( white.xyY )
             
        ACEStoObserverCAT  = spacesXYZ::CAT( white.XYZ[1, ], white.XYZ[2, ] )    # CAT from ACES whitepoint to assumed observer adapted whitepoint
        ObservertoACESCAT  = spacesXYZ::CAT( white.XYZ[2, ], white.XYZ[1, ] )    # CAT from assumed observer adapted whitepoint to ACES whitepoint
        }        
    
    surround.full   = c( 'dim', 'dark', 'normal' )
    idx = pmatch( tolower(surround), surround.full )
    if( is.na(idx) )
        {
        log_string( ERROR, "surround='%s' is invalid.", surround )
        return(NULL)
        }
    surround    = surround.full[idx]
    
    dynrange.full   = c( 'SDR', 'HDR' )
    idx = pmatch( toupper(dynrange), dynrange.full )
    if( is.na(idx) )
        {
        log_string( ERROR, "dynrange='%s' is invalid.", dynrange )
        return(NULL)
        }
    dynrange    = dynrange.full[idx]

    #   precompute a few things that are used by the 2 functions
    PARAMS  = init_TsParams3( Ymin, Ymid, Ymax )     #; print( str(PARAMS) )
    if( is.null(PARAMS) )   return(NULL)
    
    theAffine   = affine.TF( Ymin, Ymax )
        
    dataXYZ = calculateDataXYZ( display_pri, 1 )
    if( is.null(dataXYZ) )  return(NULL)
    
    XYZ_2_DISPLAY_PRI_MAT = dataXYZ$XYZ2RGB
    DISPLAY_PRI_2_XYZ_MAT = dataXYZ$RGB2XYZ
    
    if( ! is.null(limiting_pri)  &&  all( limiting_pri == display_pri ) )
        #   ignore limiting_pri, if it is the same as display_pri
        limiting_pri = NULL
        

    if( ! is.null(limiting_pri) )
        {
        datatemp = calculateDataXYZ( limiting_pri, 1 )        
        if( is.null(datatemp) )  return(NULL)        
        XYZ_2_LIMITING_PRI_MAT  = datatemp$XYZ2RGB 
        LIMITING_PRI_2_XYZ_MAT  = datatemp$RGB2XYZ 
        }        
        
        
        
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
        
        

    fun <- function( aces )
        {
        if( glowmodifier )
            #   --- Glow module ---     
            aces = glow( aces )
        
        if( redmodifier )
            #--- Red modifier ---#    
            aces = redmod( aces )
        
            
        #   if( any( aces < 0 ) )    print( aces )

        #---   AP0 to AP1  (AP1 is RGB rendering space) both are scene-linear ---# 
        #aces    = pmax( aces, 0 )      #   avoids saturated negative colors from becoming positive in the matrix
        
        rgbPost = tcrossprod( aces, p.AP0_2_AP1_MAT )
        
        #   if( any( rgbPost < 0 ) )    print( rgbPost )
        
        rgbPost = pmax( rgbPost, 0 )        # copied from rrt_sweeteners() in CTL
        
        #--- Global desaturation ---#
        rgbPost = tcrossprod( rgbPost, RRT_SAT_MAT )

        #   Apply the tonescale independently in rendering-space RGB            
        for( k in 1:3 )
            rgbPost[k] = ssts( rgbPost[k], PARAMS )     # input values are clipped to HALF_MIN, which is 5.96046448e-08
        

                
        
        #   At this point data encoded AP1, scaled absolute luminance (cd/m^2)

        #   Scale absolute luminance to linear code value 
        linearCV  = transfer( theAffine^-1, rgbPost, domaincheck=FALSE )
        
        #   Rendering primaries to XYZ
        XYZ = tcrossprod( linearCV, p.AP1_2_XYZ_MAT )

        #   Apply gamma adjustment to compensate for dim surround

        #   NOTE: This is more or less a placeholder block and is largely inactive 
        #   in its current form. This section currently only applies for SDR, and
        #   even then, only in very specific cases.
        #   In the future it is fully intended for this module to be updated to 
        #   support surround compensation regardless of luminance dynamic range. 

        #   TOD0: Come up with new surround compensation algorithm, applicable 
        #   across all dynamic ranges and supporting dark/dim/normal surround.  

        if( surround=='dim'  &&  dynrange=='SDR' )
            { 
            #   INACTIVE for HDR and crudely implemented for SDR (see comment below)        
            #   the data is SDR and so the SDR gamma compensation factor from v1.0 will apply. 
            #   This uses a local dark_to_dim function that is designed to take in
            #   XYZ and return XYZ rather than AP1 as is currently in the functions
            #   in 'ACESlib.ODT_Common.ctl' 
            XYZ = dark_to_dim( XYZ ) 
            }

        if( ! is.null(limiting_pri) )
            {
            #   sets of primaries are different
            XYZ = limit_to_primaries( XYZ, XYZ_2_LIMITING_PRI_MAT, LIMITING_PRI_2_XYZ_MAT )
            }
        

        if( ! is.null(ACEStoObserverCAT) )
            {
            #   adapt from ACES whitepoint to assumed observer adapted whitepoint
            XYZ = spacesXYZ::adaptXYZ( ACEStoObserverCAT, XYZ )   #, D60_2_D65_CAT);
            } 
        
        #   CIE XYZ to display encoding primaries, i.e. to linear display RGB
        linearCV    = tcrossprod( XYZ, XYZ_2_DISPLAY_PRI_MAT )

        #   Scale to avoid clipping when device calibration is different from D60. 
        #   To simulate D60, unequal code values are sent to the display.
        #   TODO: Needs to expand from just supporting D60 sim to allow for any
        #   observer adapted white point.
        if( is.null(ACEStoObserverCAT) )
            {
            #   observerWP is ACES D60, or unknown            
            #   TODO: The scale requires calling itself. Scale is same no matter the luminance.
            #   Currently precalculated for D65, DCI. If DCI, roll_white_fwd is used also.
            #   This needs a more complex algorithm to handle all cases.
            
            if( all( display_pri[4, ] == REC709_PRI[4, ] ) )
                {
                #   display has D65 whitepoint
                SCALE       = 0.96362
                linearCV    = SCALE * linearCV                
                }
            else if( all( display_pri[4, ] == P3DCI_PRI[4, ] ) )
                {
                #   display has DCI whitepoint
                SCALE       = 0.96                 
                linearCV    = SCALE * roll_white_fwd( linearCV, 0.918, 0.5 )
                } 
            }

        linearCV = pmax( linearCV, 0 )
        
        return( linearCV )      # now linear display RGB
        } #   end of fun()


        
        

    funinv <- function( linearCV )  #   linearCV is display linear RGB
        {
        #   Un-scale
        if( is.null(ObservertoACESCAT)  ) 
            {
            #   observerWP is ACES D60, or unknown
            #   TODO: The scale requires calling itself. Need an algorithm for this.
            #   Scale is same no matter the luminance.
            #   Currently using precalculated values for D65 and DCI.
            #   If DCI, roll_white_rev is used also.

            if( all( display_pri[4, ] == REC709_PRI[4, ] ) ) 
                {
                #   display has D65 whitepoint
                SCALE       = 0.96362
                linearCV    = linearCV / SCALE
                } 
            else if( all( display_pri[4, ] == P3DCI_PRI[4, ] ) )
                {
                #   display has DCI whitepoint
                SCALE       = 0.96
                linearCV    = roll_white_rev( linearCV/SCALE, 0.918, 0.5)
                }
            }
        

        #   Encoding primaries to CIE XYZ
        XYZ = tcrossprod( linearCV, DISPLAY_PRI_2_XYZ_MAT )


        if( ! is.null(ObservertoACESCAT) )
            {
            #   print( XYZ ) ; print( dim(XYZ) )            
            XYZ = spacesXYZ::adaptXYZ( ObservertoACESCAT, XYZ )    # adapt from display (assumed observer adapted) whitepoint to ACES whitepoint
            }
            
        if( ! is.null(limiting_pri) )
            {
            #   primaries are different
            #   what goes here ?   How do you unlimit ?
            }            
            
        #   unapply gamma adjustment to compensate for dim surround

        #   NOTE: This is more or less a placeholder block and is largely inactive 
        #   in its current form. This section currently only applies for SDR, and
        #   even then, only in very specific cases.
        #   In the future it is fully intended for this module to be updated to 
        #   support surround compensation regardless of luminance dynamic range. 

        #   TOD0: Come up with new surround compensation algorithm, applicable 
        #   across all dynamic ranges and supporting dark/dim/normal surround.  

        if( surround=='dim'  &&  dynrange=='SDR' )
            { 
            #   INACTIVE for HDR and crudely implemented for SDR (see comment below)        
            #   the data is SDR and so the SDR gamma compensation factor from v1.0 will apply. 
            #   This uses a local dark_to_dim function that is designed to take in
            #   XYZ and return XYZ rather than AP1 as is currently in the functions
            #   in 'ACESlib.ODT_Common.ctl' 
            XYZ = dim_to_dark( XYZ ) 
            }
            
        #   XYZ to rendering primaries
        linearCV = tcrossprod( XYZ, p.XYZ_2_AP1_MAT )
            
        #if( any( 1<linearCV ) )    { print( linearCV ) }
            
        #linearCV = pmax( linearCV, 0 )    #; print( linearCV )  # not in original aces-dev

        #   Scale linear code value to absolute luminance
        rgbPost = transfer( theAffine, linearCV, domaincheck=FALSE )  #;        print( rgbPost )
        
        #   Apply the inverse tonescale independently in rendering-space RGB             
        #   rgbPre  = rgbPost
        for( k in 1:3 )
            rgbPost[k]   = inv_ssts( rgbPost[k], PARAMS )

        #   print( rgbPost )
        
        #--- Global desaturation ---#
        rgbPost = tcrossprod( rgbPost, RRT_SAT_MAT_inv )

        #--- AP1 (ACES RGB rendering space) to AP0 ---#    
        #   rgbPost = pmax( rgbPost, 0 )    #clamp_f3( rgbPost, 0., HALF_MAX);
        aces    = tcrossprod( rgbPost, p.AP1_2_AP0_MAT )
        #   aces    = rgbPost
        #   aces    = pmax( aces, 0 )      #   avoids saturated negative colors from becoming positive in the matrix
        
        if( redmodifier )
            {
            #--- Red modifier inverse   AP0 ---#
            if( preciseinv )
                aces = redmodinv_precise( aces )
            else
                aces = redmodinv( aces )    # rough approximation
            }
        
        if( glowmodifier )        
            #--- Glow module inverse ---#
            aces = glowinv( aces )
    
        return( aces )
        }   #   end of funinv()
    
    rgbinterval = c( PARAMS$Min$x, PARAMS$Max$x )   #1 )
    
    cnames  = sprintf( "ACES.%s", c('R','G','B') )
    domain  = matrix( rgbinterval, 2, 3, dimnames=list(NULL,cnames) )
    
    #   make all vertices of cube, in 8x3 matrix
    mat     = as.matrix( expand.grid( rgbinterval, rgbinterval, rgbinterval ) )
    for( i in 1:nrow(mat) )
        mat[i, ] = fun( mat[i, ] )
    
    orange              = apply( mat, 2, range ) #; print(orange)
    orange              = matrix( c(0,1), 2, 3 )
    colnames(orange)    = sprintf( "displaylinear.%s", c('R','G','B') )    

    out = TransferFunction( fun, funinv, domain, orange, id=sigfunction() )
    
    #   in the next line, the variable white=Ymax may be used in installRGB
    metadata(out)   = list( primaries=dataXYZ$primaries, white=Ymax )   # dataXYZ$whiteXYZ[2] )
    
    return( out )
    }   #   end of general.OOTF()
    
    
    
dark_to_dim <-  function( XYZ )
    {
    if( XYZ[2] <= 0 )   return( numeric(3) )
    
    s   = XYZ[2]^(DIM_SURROUND_GAMMA - 1)
    
    return( s * XYZ )
    }

dim_to_dark <-  function( XYZ )
    {
    if( XYZ[2] <= 0 )   return( numeric(3) )
    
    s   = XYZ[2]^(1/DIM_SURROUND_GAMMA - 1)

    return( s * XYZ )
    }
    
    
#   XYZ  ->  RGB
#   clamp RGB to unit cube
#   RGB  ->  XYZ

limit_to_primaries  <- function( XYZ, XYZ_2_LIMITING_PRI_MAT, LIMITING_PRI_2_XYZ_MAT )
    {
    #   XYZ to limiting primaries
    RGB = tcrossprod( XYZ, XYZ_2_LIMITING_PRI_MAT )

    #   Clip any values outside the limiting primaries
    limitedRgb  = pmin( pmax(RGB,0), 1 )
    
    #   Convert limited RGB to XYZ
    XYZ = tcrossprod( limitedRgb, LIMITING_PRI_2_XYZ_MAT )
    
    return( XYZ )
    }
    