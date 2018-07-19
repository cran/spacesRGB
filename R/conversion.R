

###################################      XYZ  <-->  RGB    ################################

#   RGB         n x 3 matrix
#   space       RGB space, e.g. 'sRGB', 'AdobeRGB', ...
#   which       which linear RGB, 'scene' or 'display'
#   TF          transfer function, or numeric gamma  
#   maxSignal   1, 255, 1023, etc.
#
#   returns a data.frame with 2 columns
#           XYZ         either scene XYZ or display XYZ
#           OutOfGamut  logical
#   the returned XYZ is for viewing under the white-point of the RGB space (usually D65)


XYZfromRGB <- function( RGB, space='sRGB', which='scene', TF=NULL, maxSignal=1 )
    {
    out = LinearRGBfromSignalRGB( RGB, space=space, which=which, TF=TF, maxSignal=maxSignal )
    
    if( is.null(out) )  return(NULL)
                
    # verify space                          
    idx = spaceIndex( space )
    if( is.na(idx) )    return(NULL)

    out$XYZ = tcrossprod( out$RGB, p.ListRGB[[idx]]$RGB2XYZ )
    #   colnames(out$XYZ) = c('X','Y','Z')
    
    out$RGB = NULL
    
    return( out )
    }


#   XYZ         for viewing under the white-point of the RGB space (usually D65), with Y=100 for white
#   space       RGB space, e.g. 'sRGB', 'AdobeRGB', ...
#   which       which linear RGB, 'scene' or 'display'
#   TF          transfer function, or numeric gamma 
#   maxSignal   of SignalRGB
#
#   returns a data.frame with 2 columns
#           RGB         signal RGB for the device.  clamped to apropriate cube, unless gamma=1
#           OutOfGamut  logical
RGBfromXYZ <- function( XYZ, space='sRGB', which='scene', TF=NULL, maxSignal=1 )
    {
    # verify XYZ
    XYZ = prepareNxM(XYZ)
    if( is.null(XYZ) )  return(NULL)

    # verify space 
    idx = spaceIndex( space )
    if( is.na(idx) )    return(NULL)

    RGB = tcrossprod( XYZ, p.ListRGB[[idx]]$XYZ2RGB )   # t(M %*% t(XYZ))    # print( RGB[1, ] - 1 )

    return( SignalRGBfromLinearRGB( RGB, space=space, which=which, TF=TF, maxSignal=maxSignal ) )
    }



#   RGB         linear RGB any sort of array
#               it is OK if val is a matrix, and then the return value is a matrix of the same shape
#   space       name of color space
#   which       which linear RGB, 'scene' or 'display'
#   TF          transfer function, or numeric gamma of the display, so output is (linear)^(1/gamma)
#   maxSignal   of SignalRGB
#
#   return  first clips to [0,1], and then maps [0,1] to [0,1].
#           in case of ERROR it logs a message and returns the clipped values only
#
#   returns a data.frame with 2 columns
#           RGB         signal RGB for the device.  clamped to apropriate cube, unless gamma=1
#           OutOfGamut  logical, TRUE if clamping actually performed

SignalRGBfromLinearRGB <- function( RGB, space='sRGB', which='scene', TF=NULL, maxSignal=1 )
    {
    # verify RGB
    RGB = prepareNxM(RGB)
    if( is.null(RGB) )  return(NULL)

    # The out-of-gamut flag is a column vector of Boolean true/false values.  Each
    # entry corresponds to one row of the input matrix XYZ.
    # There is numerical tolerance here, designed for points in XYZ that should map to vertices of the RGB cube
    lo  = -1.e-7
    hi  = 1 + 1.e-7
    OutOfGamut  = (RGB[ ,1] < lo  |  RGB[ ,1] > hi  |  RGB[ ,2] < lo  |  RGB[ ,2] > hi  |  RGB[ ,3] < lo  |  RGB[ ,3] > hi)

    
    if( is.null(TF) )
        {
        # verify space 
        idx = spaceIndex( space )
        if( is.na(idx) )    return(NULL)
            
        theSpace   = p.ListRGB[[idx]]
    
        if( which == 'scene' )
            TF  = theSpace$OETF         # forward to signal
        else if( which == 'display' )
            TF  = theSpace$EOTFinv      # backward to signal
        else
            {
            log.string( ERROR, "which='%s' is invalid.", which )
            return(NULL)
            }
        }

    # convert from linear to signal
    if( is.function(TF) )
        {
        #   clamp RGB to the unit cube  
        RGB[ RGB<0 ] = 0
        RGB[ RGB>1 ] = 1 
        
        RGBsignal  =  TF( RGB )    
        }
    else if( is.numeric(TF)  &&  length(TF)==1  &&  0 < TF )
        {
        if( TF == 1 )
            RGBsignal   = RGB   # no clamping, identity transfer
        else
            {
            #   clamp RGB to the unit cube  
            RGB[ RGB<0 ] = 0
            RGB[ RGB>1 ] = 1
            RGBsignal    = RGB ^ (1/TF)
            }
        }
    else
        {
        log.string( ERROR, "TF is invalid."  )
        return(NULL)
        }

    RGBsignal = maxSignal * RGBsignal

    colnames(RGBsignal) = c('R','G','B')

    rnames  = rownames(RGB)
    if( is.null(rnames)  ||  0<anyDuplicated(rnames)  )
        #   rnames is no good !  Use trivial names instead.
        rnames = 1:nrow(RGB)

    out = data.frame( row.names=rnames )
    out$RGB             = RGBsignal
    out$OutOfGamut      = OutOfGamut 

    return( out )
    }
    
    
#   returns a data.frame with 2 columns
#           RGB         linear RGB, either 'scene' or 'display'.  clamped to apropriate cube, unless gamma=1
#           OutOfGamut  logical, TRUE if clamping actually performed
    
LinearRGBfromSignalRGB <- function( RGB, space='sRGB', which='scene', TF=NULL, maxSignal=1 )
    {
    # verify RGB
    RGB = prepareNxM(RGB)
    if( is.null(RGB) )  return(NULL)
    
    # verify maxSignal
    ok  = is.numeric(maxSignal)  &&  length(maxSignal)==1  &&  0<maxSignal
    if( ! ok )
        {
        log.string( ERROR, "maxSignal='%s' is not a positive number.", as.character(maxSignal) )
        return(NULL)
        }    
                         
    RGB = RGB/maxSignal
    
    lo  = 0
    hi  = 1
    OutOfGamut  = (RGB[ ,1] < lo  |  RGB[ ,1] > hi  |  RGB[ ,2] < lo  |  RGB[ ,2] > hi  |  RGB[ ,3] < lo  |  RGB[ ,3] > hi)
    
    if( is.null(TF) )
        {
        # verify space 
        idx = spaceIndex( space )
        if( is.na(idx) )    return(NULL)

        theSpace   = p.ListRGB[[idx]]

        if( which == 'scene' )
            TF  = theSpace$OETFinv      # backwards to scene
        else if( which == 'display' )
            TF  = theSpace$EOTF         # forwards to display
        else
            {
            log.string( ERROR, "which='%s' is invalid.", as.character(which) )
            return(NULL)
            }            
        }

    if( is.function( TF ) )
        {
        #   clamp RGB to the unit cube (inside gamut)
        RGB[ RGB<0 ]    = 0
        RGB[ RGB>1 ]    = 1         
        
        RGBlin  = TF( RGB ) 
        }
    else if( is.numeric(TF) &&  length(TF)==1  &&  0 < TF )
        {
        if( TF == 1 )
            RGBlin = RGB        # identity, so copy without clamping
        else
            {
            #   clamp RGB to the unit cube (inside gamut)
            RGB[ RGB<0 ]    = 0
            RGB[ RGB>1 ]    = 1
            
            RGBlin  = RGB ^ TF
            }
        }
    else
        {
        log.string( ERROR, "argument TF is invalid."  )
        return(NULL)
        }
    
    colnames(RGBlin) = c('R','G','B')
    
    rnames  = rownames(RGB)
    if( is.null(rnames)  ||  0<anyDuplicated(rnames)  )
        #   rnames is no good !  Use trivial names instead.
        rnames = 1:nrow(RGB)

    out = data.frame( row.names=rnames )
    out$RGB             = RGBlin
    out$OutOfGamut      = OutOfGamut
    
    return( out )
    }
    


