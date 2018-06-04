

###################################      XYZ  <-->  RGB    ################################

#   RGB         n x 3 matrix
#   space       RGB space, e.g. 'sRGB', 'AdobeRGB', ...
#   gamma       if not NULL, custom gamma to override the gamma of space.  If gamma is 1, input RGB is linear
#   maxValue    1, 255, 1023, etc.
#
#   the returned XYZ is for viewing under the white-point of the RGB space (usually D65)
#   and when RGB are all maxValue, then Y=100

XYZfromRGB <- function( RGB, space='sRGB', gamma=NULL, maxValue=1 )
    {
    out = LinearRGBfromDisplayRGB( RGB, space=space, gamma=gamma, maxValue=maxValue )
    
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
#   gamma       if not NULL, custom gamma to override the gamma of space.  If gamma is 1, output RGB is linear.
#   maxValue    for the device.  1023 is also common
#
#   returns a data.frame with 2 columns
#           RGB             for the device.  clamped to apropriate cube, unless gamma=1
#           OutOfGamutFlag  logical
RGBfromXYZ <- function( XYZ, space='sRGB', gamma=NULL, maxValue=1 )
    {
    # verify XYZ
    XYZ = prepareNxM(XYZ)
    if( is.null(XYZ) )  return(NULL)

    # verify space 
    idx = spaceIndex( space )
    if( is.na(idx) )    return(NULL)

    RGB = tcrossprod( XYZ, p.ListRGB[[idx]]$XYZ2RGB )   # t(M %*% t(XYZ))    # print( RGB[1, ] - 1 )

    return( DisplayRGBfromLinearRGB( RGB, space=space, gamma=gamma, maxValue=maxValue ) )
    }



#   RGB     linear RGB any sort of array
#           it is OK if val is a matrix, and then the return value is a matrix of the same shape
#   space   name of color space
#
#   gamma   transfer function, or numeric gamma of the display, so output is (linear)^(1/gamma)
#
#   maxValue    of DisplayRGB
#
#   return  first clips to [0,1], and then maps [0,1] to [0,1].
#           in case of ERROR it logs a message and returns the clipped values only
#
DisplayRGBfromLinearRGB <- function( RGB, space='sRGB', gamma=NULL, maxValue=1 )
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

    
    if( is.null(gamma) )
        {
        # verify space 
        idx = spaceIndex( space )
        if( is.na(idx) )    return(NULL)
            
        theSpace   = p.ListRGB[[idx]]
    
        if( is.list(theSpace$gamma) )
            gamma = theSpace$gamma[[2]]     # a function
        else
            gamma = theSpace$gamma          # a number
        }

    # convert from linear to display
    if( is.function( gamma ) )
        {
        #   clamp RGB to the unit cube  
        RGB[ RGB<0 ] = 0
        RGB[ RGB>1 ] = 1 
        
        RGBdisplay  =  gamma( RGB )    
        }
    else if( is.numeric(gamma) &&  length(gamma)==1  &&  0 < gamma )
        {
        if( gamma == 1 )
            RGBdisplay  = RGB   # no clamping
        else
            {
            #   clamp RGB to the unit cube  
            RGB[ RGB<0 ] <- 0
            RGB[ RGB>1 ] <- 1
            
            RGBdisplay  = RGB ^ (1/gamma)
            }
        }
    else
        {
        log.string( ERROR, "gamma is invalid."  )
        return(NULL)
        }

    RGBdisplay = maxValue * RGBdisplay

    colnames(RGBdisplay) = c('R','G','B')

    rnames  = rownames(RGB)
    if( is.null(rnames)  ||  0<anyDuplicated(rnames)  )
        #   rnames is no good !  Use trivial names instead.
        rnames = 1:nrow(RGB)

    out = data.frame( row.names=rnames )
    out$RGB             = RGBdisplay
    out$OutOfGamut      = OutOfGamut 

    return( out )
    }
    
    
LinearRGBfromDisplayRGB <- function( RGB, space='sRGB', gamma=NULL, maxValue=1 )
    {
    # verify RGB
    RGB = prepareNxM(RGB)
    if( is.null(RGB) )  return(NULL)
    
    # verify maxValue
    ok  = is.numeric(maxValue)  &&  length(maxValue)==1  &&  0<maxValue
    if( ! ok )
        {
        log.string( ERROR, "maxValue='%s' is not a positive number.", as.character(maxValue) )
        return(NULL)
        }    
                         
    RGB = RGB/maxValue
    
    lo  = 0
    hi  = 1
    OutOfGamut  = (RGB[ ,1] < lo  |  RGB[ ,1] > hi  |  RGB[ ,2] < lo  |  RGB[ ,2] > hi  |  RGB[ ,3] < lo  |  RGB[ ,3] > hi)
    
    if( is.null(gamma) )
        {
        # verify space 
        idx = spaceIndex( space )
        if( is.na(idx) )    return(NULL)

        theSpace   = p.ListRGB[[idx]]

        if( is.list(theSpace$gamma) )
            gamma = theSpace$gamma[[1]]     # a function
        else
            gamma = theSpace$gamma          # a number
        }

    if( is.function( gamma ) )
        {
        #   clamp RGB to the unit cube (inside gamut)
        RGB[ RGB<0 ]    = 0
        RGB[ RGB>1 ]    = 1         
        
        RGBlin  = gamma( RGB ) 
        }
    else if( is.numeric(gamma) &&  length(gamma)==1  &&  0 < gamma )
        {
        if( gamma == 1 )
            RGBlin = RGB
        else
            {
            #   clamp RGB to the unit cube (inside gamut)
            RGB[ RGB<0 ]    = 0
            RGB[ RGB>1 ]    = 1
            
            RGBlin  = RGB ^ gamma
            }
        }
    else
        {
        log.string( ERROR, "argument gamma is invalid."  )
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
    

#   space           name of the space
#   return value    index in p.ListRGB, partial matching and case-insensitive
#                   if no match, then NA_integer_
spaceIndex <- function( space )
    {
    ok  = is.character(space)  &&  length(space)==1
    if( ! ok )
        {
        log.string( ERROR, "space is not a character vector of length 1." )
        return(NA_integer_)
        }
        
    theNames    = names(p.ListRGB)
    if( is.null(theNames)  ||  length(theNames)==0 )
        {
        log.string( ERROR, "ERROR internal.  There are no installed RGB spaces." )
        return(NA_integer_)
        }

    idx = pmatch( toupper(space), toupper(theNames) )
    if( is.na(idx) )
        {
        log.string( ERROR, "space='%s' matches no installed spaces, or multiple spaces.", space )
        return(NA_integer_)
        }

    return(idx)
    }
    
