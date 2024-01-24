
#   PODT =  Partial ODT = Partial Output Device Transform
#
#   general.PODT() is missing the final OETF; the output is still 'linear' 


#   Target white and black points for cinema system tonescale
CINEMA_WHITE    = 48
CINEMA_BLACK    = 10 ^ log10(0.02)  # CINEMA_WHITE / 2400.   ~=~ 0.02 +- 1.e-18
#   CINEMA_BLACK is defined in this roundabout manner in order to be exactly equal to 
#   the result returned by the cinema 48-nit ODT tonescale.
#   Though the min point of the tonescale is designed to return 0.02, the tonescale is 
#   applied in log-log space, which loses precision on the antilog. The tonescale 
#   return value is passed into Y_2_linCV, where CINEMA_BLACK is subtracted. If 
#   CINEMA_BLACK is defined as simply 0.02, then the return value of this subfunction
#   is very, very small but not equal to 0, and attaining a CV of 0 is then impossible.
#   For all intents and purposes, CINEMA_BLACK=0.02.


#   display_pri     a 4x2 matrix with xy chromaticities of RGBW in the rows.  The display primaries.
#                   If this is NULL, as in the case of DCDM, it means that the output RGB is really XYZ
#   Ymax            maximum luminance, in nits.  This is only stored in the metadata of the returned TransferFunction.
#   observerWP      assumed observer adapted-whitepoint.  This is used to form a CAT that goes from D60 to the observer WP.
#                   If this is NULL, then the function will get it from display_pri.
#                   If display_pri is NULL, then the function will get it from limiting_pri.
#                   And limiting_pri is also NULL, there is no CAT.
#                   And if the observer WP is D60, there is no CAT.
#   surround        'dark', or 'dim'. If 'dark' then do nothing.  If 'dim' then compensate and also desaturate.
#   limiting_pri    a 4x2 matrix with xy chromaticities of RGBW in the rows
#                   If this is NULL, there is no limiting (clipping) to the unit RGB cube (limiting RGB).

general.PODT <- function(   display_pri, Ymax=1,
                            observerWP=NULL,
                            surround='dark',
                            limiting_pri=NULL )
    {
    if( ! is.null(display_pri) )
        {
        dataXYZ = calculateDataXYZ( display_pri, 1 )
        if( is.null(dataXYZ) )  return(NULL)
        
        XYZ_2_DISPLAY_PRI_MAT = dataXYZ$XYZ2RGB
        DISPLAY_PRI_2_XYZ_MAT = dataXYZ$RGB2XYZ
        }
        
    if( is.null(observerWP) )
        {
        #   try to get observerWP from the primaries
        if( ! is.null(display_pri) )
            observerWP  = display_pri[4, ]
        else if( ! is.null(limiting_pri) )
            observerWP  = limiting_pri[4, ]
        }
        
        
    if( is.null(observerWP)   ||  all( observerWP == AP0_PRI[4, ] ) )       # ||  ! is.na(simscale) )
        {
        #   observer whitepoint is not available, or it matches AP0 whitepoint  - so no need to make CAT
        ACEStoObserverCAT = NULL
        ObservertoACESCAT = NULL
        }    
    else
        {
        #   TODO: Needs to expand from just supporting D60 sim to allow for any observer-adapted white point.
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

    if( ! is.null(limiting_pri) )
        {
        datatemp = calculateDataXYZ( limiting_pri, 1 )        
        XYZ_2_LIMITING_PRI_MAT  = datatemp$XYZ2RGB 
        LIMITING_PRI_2_XYZ_MAT  = datatemp$RGB2XYZ 
        }
        
    #   process surround argument
    surround.full   = c( 'dim', 'dark' )    #, 'normal' )
    idx = pmatch( tolower(surround), surround.full )
    if( is.na(idx) )
        {        
        log_string( ERROR, "surround='%s' is invalid.", surround )
        return(NULL)
        }
    surround    = surround.full[idx]        
        
    #   assign simscale and rolloff from display_pri and observerWP
    rolloff     = rep( NA_real_, 2 )                           
    simscale    = NA_real_
        
    if( identical(display_pri,P3D65_PRI)  &&  identical(observerWP,AP0_PRI['W', ]) )
        {
        simscale    = 0.964
        }
    else if( identical(display_pri,P3DCI_PRI) )
        {
        if( identical(observerWP,AP0_PRI['W', ]) )
            {
            simscale    = 0.96
            rolloff     = c( 0.918, 0.5 )
            }
        else if( identical(observerWP,P3D65_PRI['W', ]) )
            {
            simscale    = 0.9575
            rolloff     = c( 0.908, 0.5 )
            }
        }
    else if( identical(display_pri,REC709_PRI)  &&  identical(observerWP,AP0_PRI['W', ]) )
        {
        simscale    = 0.955
        }
        
    #print( simscale )
    #print( rolloff )
        
        
    fun <- function( oces )
        {
        #   OCES to RGB rendering space - AP1
        rgbPre  =  tcrossprod( oces, p.AP0_2_AP1_MAT )   #; cat( 'rgbPre=', rgbPre, '\n' )

        #   Apply the tonescale independently in rendering-space RGB
        rgbPost = base::sapply( rgbPre, segmented_spline_fwd, C=p.ODT_48nits )  #; cat( 'rgbPost=', rgbPost, '\n' )

        
        #   Scale luminance to linear code value
        linearCV    = (rgbPost - CINEMA_BLACK) / (CINEMA_WHITE - CINEMA_BLACK)  #; cat( 'linearCV=', linearCV, '\n' )

        if( all( ! is.na(rolloff) ) )
            {
            #   However, the magnitude of the scale factor required was considered too 
            #   large; therefore, the scale factor was reduced and the additional 
            #   required compression was achieved via a reshaping of the highlight 
            #   rolloff in conjunction with the scale. The shape of this rolloff was 
            #   determined through subjective experiments and deemed to best reproduce 
            #   the "character" of the highlights in the P3D60 ODT.
    
            #   Roll off highlights to avoid need for as much scaling
            NEW_WHT     = rolloff[1]
            ROLL_WIDTH  = rolloff[2]
            linearCV    = roll_white_fwd( linearCV, NEW_WHT, ROLL_WIDTH )
            }
        
        if( ! is.na(simscale)  )
            {
            #   --- Compensate for different white point being darker  --- //
            #   This adjustment corrects for an issue that exists in ODTs where the 
            #   device is calibrated to a white chromaticity other than that of the adapted white.
            #   In order to produce D60 on a device calibrated to D65 white point (i.e. 
            #   equal display code values yield CIE x,y chromaticities of 0.3127, 0.329) 
            #   the red channel is higher than green and blue to compensate for the 
            #   "bluer" D65 white. This is the intended behavior but it means that 
            #   without compensation, as highlights increase, the red channel will hit 
            #   the device maximum first and clip, resulting in a chromaticity shift as 
            #   the green and blue channels continue to increase.
            #   To avoid this clipping behavior, a slight scale factor is applied to 
            #   allow the ODTs to simulate D60 within the D65 calibration white point. 

            #   scale and clamp white to avoid casted highlights due to D60 simulation
            linearCV    = pmin( linearCV, 1 ) * simscale       
            }

            
        if( surround == 'dim' )
            {
            #   Apply gamma adjustment to compensate for dim surround            
            #   TOD0: Come up with new surround compensation algorithm, applicable 
            #   across all dynamic ranges and supporting dark/dim/normal surround.              
            linearCV    = darkSurround_to_dimSurround( linearCV )   # this goes to XYZ, modifies XYZ, and then back again !

            #   Apply desaturation to compensate for luminance difference
            linearCV    = tcrossprod( linearCV, ODT_SAT_MAT )
            }
            
            
        #   Convert to display primary encoding rendering space RGB - AP1,  to XYZ
        XYZ = tcrossprod( linearCV, p.AP1_2_XYZ_MAT )

        
        if( ! is.null(ACEStoObserverCAT) ) #! D60_sim  &&  
            {
            #   adapt from ACES whitepoint to assumed observer adapted whitepoint
            #   cat( "before CAT, XYZ = ", nicevector(XYZ), "   xyY = ", nicevector(xyYfromXYZ(XYZ)), '\n' )
            XYZ = spacesXYZ::adaptXYZ( ACEStoObserverCAT, XYZ )          #; print( 'adapted' )
            #   cat( "after CAT,  XYZ = ", nicevector(XYZ), "   xyY = ", nicevector(xyYfromXYZ(XYZ)), '\n' )            
            } 

        if( ! is.null(limiting_pri) )
            {
            #   clamp to limiting primaries
            XYZ = limit_to_primaries( XYZ, XYZ_2_LIMITING_PRI_MAT, LIMITING_PRI_2_XYZ_MAT ) #; print( 'limited' )
            }            
            
            
        if( ! is.null(display_pri) )
            #   CIE XYZ to display primaries
            linearCV = tcrossprod( XYZ, XYZ_2_DISPLAY_PRI_MAT )
        else
            linearCV = XYZ
            

        #   Handle out-of-gamut values
        #   Clip values < 0 or > 1 (i.e. projecting outside the display primaries)
        linearCV    = pmin( pmax(linearCV,0), 1 )   

        
        return( linearCV )
        }

    funinv <- function( linearCV )
        {        
        if( ! is.null(display_pri) )        
            #   Convert from display primary encoding display primaries to CIE XYZ
            XYZ = tcrossprod( linearCV, DISPLAY_PRI_2_XYZ_MAT )
        else
            XYZ = linearCV
            
            
        #   unlimit the limiting primaries would go here, but one cannot unlimit.
            
        if( ! is.null(ObservertoACESCAT) )   #! D60_sim  && 
            {
            #   print( XYZ ) ; print( dim(XYZ) )            
            XYZ = spacesXYZ::adaptXYZ( ObservertoACESCAT, XYZ )    # adapt from assumed observer adapted whitepoint to ACES whitepoint
            }        

        #   CIE XYZ to rendering space RGB = AP1
        linearCV    = tcrossprod( XYZ, p.XYZ_2_AP1_MAT )
        
        if( surround == 'dim' )
            {
            #   Undo desaturation to compensate for luminance difference
            linearCV    = tcrossprod( linearCV, ODT_SAT_MAT_inv )

            #   Undo gamma adjustment to compensate for dim surround
            linearCV    = dimSurround_to_darkSurround( linearCV )
            }
        
        
        
        if( ! is.na(simscale) )
            {
            #   Undo scaling done for D60 simulation
            linearCV  = linearCV  / simscale
            }
        
        if( all( ! is.na(rolloff) ) )
            {        
            #   Undo highlight roll-off and scaling
            NEW_WHT     = rolloff[1]
            ROLL_WIDTH  = rolloff[2]            
            linearCV    = roll_white_rev( linearCV, NEW_WHT, ROLL_WIDTH )
            }
            
        #   Scale linear code value to luminance
        rgbPre  = (1 - linearCV)*CINEMA_BLACK  +  linearCV*CINEMA_WHITE

        #   Apply the tonescale independently in rendering-space RGB = AP1
        rgbPost = base::sapply( rgbPre, segmented_spline_rev, C=p.ODT_48nits )
        
        #   Rendering space RGB (AP1) to OCES
        oces    = tcrossprod( rgbPost, p.AP1_2_AP0_MAT )

        return( oces )
        }
        
    if( ! is.null(limiting_pri) )
        #   override funinv, so the TransferFunction is NOT invertible
        funinv  = NULL

        
    domain  = matrix( c(0,10000), 2, 3, dimnames=list(NULL,c('OCES.R','OCES.G','OCES.B')) )
    range   = matrix( c(0,1), 2, 3, dimnames=list(NULL,c('signal.R','signal.G','signal.B')) )    

    out = TransferFunction( fun, funinv, domain, range, id=sigfunction() )        
    
    metadata(out)   = list( primaries=display_pri, white=Ymax )
        
    return( out )
    } 
    
    
#   linearCV    in AP1
darkSurround_to_dimSurround <- function( linearCV )
    {
    XYZ = tcrossprod( linearCV, p.AP1_2_XYZ_MAT )

    XYZ = dark_to_dim( XYZ )

    return( as.numeric( tcrossprod( XYZ, p.XYZ_2_AP1_MAT ) ) )
    }
    
    
#   linearCV    in AP1
dimSurround_to_darkSurround <- function( linearCV )
    {
    XYZ = tcrossprod( linearCV, p.AP1_2_XYZ_MAT )

    XYZ = dim_to_dark( XYZ )

    return( as.numeric( tcrossprod( XYZ, p.XYZ_2_AP1_MAT) ) )
    }
    