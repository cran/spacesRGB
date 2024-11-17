

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
    
    # verify which
    w   = endIndex(which)
    if( is.na(w) )      return(NULL)

    if( w == 1 )
        RGB2XYZ = p.ListRGB[[idx]]$scene$RGB2XYZ
    else
        RGB2XYZ = p.ListRGB[[idx]]$display$RGB2XYZ
    
    out$XYZ = tcrossprod( out$RGB, RGB2XYZ )
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
    
    # verify which
    w   = endIndex(which)
    if( is.na(w) )      return(NULL)

    if( w == 1 )
        XYZ2RGB = p.ListRGB[[idx]]$scene$XYZ2RGB
    else
        XYZ2RGB = p.ListRGB[[idx]]$display$XYZ2RGB    

    RGB = tcrossprod( XYZ, XYZ2RGB )   # t(M %*% t(XYZ))    # print( RGB[1, ] - 1 )

    return( SignalRGBfromLinearRGB( RGB, space=space, which=which, TF=TF, maxSignal=maxSignal ) )
    }


#   converts RGB -> XYZ  and then  XYZ -> Lab
#   special attention is given to the case when R==G==B, and then a=b=0

LabfromRGB <- function( RGB, space='sRGB', which='scene', TF=NULL, maxSignal=1 )
    {
    if( ! requireNamespace( 'spacesXYZ', quietly=TRUE ) )
        {    
        log_string( ERROR, "'spacesXYZ' cannot be loaded." )
        return( NULL )
        }
        
    # verify RGB
    RGB = prepareNxM(RGB)
    if( is.null(RGB) )  return(NULL)        
            
    out = XYZfromRGB( RGB, space=space, which=which, TF=TF, maxSignal=maxSignal )
    if( is.null(out) )  return(NULL)
    
    white   = getWhiteXYZ( space, which=which )
    
    Lab = spacesXYZ::LabfromXYZ( out$XYZ, white )
    
    #   for exact neutrals, set a=b=0 exactly
    neutral = RGB[ ,1]==RGB[ ,2]  &  RGB[ ,2]==RGB[ ,3]
    if( any(neutral) )
        Lab[neutral,2:3]  = 0
        
    #   change the name of the 'XYZ' column to 'Lab', and overwrite it with Lab
    colnames(out)['XYZ'==colnames(out)] = 'Lab'
    out$Lab = Lab

    return( out )
    }


RGBfromLab <- function( Lab, space='sRGB', which='scene', TF=NULL, maxSignal=1 )
    {
    if( ! requireNamespace( 'spacesXYZ', quietly=TRUE ) )
        {    
        log_string( ERROR, "'spacesXYZ' cannot be loaded." )
        return( NULL )
        }
        
    # verify Lab
    Lab = prepareNxM(Lab)
    if( is.null(Lab) )  return(NULL)

    # verify space 
    idx = spaceIndex( space )
    if( is.na(idx) )    return(NULL)
    
    # verify which
    w   = endIndex(which)
    if( is.na(w) )      return(NULL)
    
    white   = getWhiteXYZ( space, which=which )
        
    XYZ = spacesXYZ::XYZfromLab( Lab, white=white )
    
    if( w == 1 )
        XYZ2RGB = p.ListRGB[[idx]]$scene$XYZ2RGB
    else
        XYZ2RGB = p.ListRGB[[idx]]$display$XYZ2RGB    

    RGB = tcrossprod( XYZ, XYZ2RGB )   # t(M %*% t(XYZ))    # print( RGB[1, ] - 1 )
    
    #   for exact neutrals, set R=G=B = the average of R,G,B 
    neutral = Lab[ ,2]==0  &  Lab[ ,3]==0
    if( any(neutral) )
        RGB[neutral, ] = rowMeans( RGB[neutral, ,drop=FALSE] )

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

    # verify which
    w   = endIndex(which)
    if( is.na(w) )      return(NULL)    
    

    if( is.null(TF) )
        {
        # verify space 
        idx = spaceIndex( space )
        if( is.na(idx) )    return(NULL)

        theSpace   = p.ListRGB[[idx]]
    
        if( w == 1 )
            TF  = theSpace$OETF         # forwards from scene to signal
        else if( is.invertible(theSpace$EOTF) )
            TF  = (theSpace$EOTF)^-1    # backwards from display to signal
        else 
            {
            log_string( ERROR, "For RGB space='%s', the EOTF is not invertible.", space )
            return(NULL)
            }    
        }
    else if( is.numeric(TF)  &&  length(TF)==1  &&  0 < TF )
        {
        if( TF == 1 )
            TF  = identity.TF           # no clamping
        else if( w == 1 )
            TF  = power.OETF(TF)        # forward from scene to signal
        else
            TF  = power.EOTF(TF)^-1     # backward from display  to signal
        }
        
    if( ! is.TransferFunction(TF) )        
        {
        log_string( ERROR, "argument TF is invalid."  )
        return(NULL)
        }        
        
    rnames  = rownames(RGB)
    if( is.null(rnames)  ||  0<anyDuplicated(rnames)  )
        #   rnames is no good !  Use trivial names instead.
        rnames = 1:nrow(RGB)

    #   get the domain of TF and set lo and hi from that domain        
    if( is.identity(TF) )
        {
        #    no domain available, use unit cube
        lo  = 0
        hi  = 1       
        tol = 1.e-7
        OutOfGamut  = (RGB[ ,1] < lo-tol  |  RGB[ ,1] > hi+tol  |  RGB[ ,2] < lo-tol  |  RGB[ ,2] > hi+tol  |  RGB[ ,3] < lo-tol  |  RGB[ ,3] > hi+tol)

        #   no clamping
        RGBsignal   = RGB

        #   OutOfGamut  = FALSE     # logical( nrow(RGB) )
        }
    else
        {
        # convert from linear to signal        
        domain  = domain(TF)    # TF$element[[1]]$domain

        lo      = domain[1, ]       #lo  = -1.e-7
        hi      = domain[2, ]       #hi  = 1 + 1.e-7        
        
        if( length(lo) == 1 )
            {
            lo = rep(lo,3)
            hi = rep(hi,3)
            }
            
        # The out-of-gamut flag is a column vector of Boolean true/false values.  Each
        # entry corresponds to one row of the input matrix RGB.
        # There is numerical tolerance here, designed for points in XYZ that should map to vertices of the RGB cube
        tol = 1.e-7
        OutOfGamut  = (RGB[ ,1] < lo[1]-tol  |  RGB[ ,1] > hi[1]+tol  |  RGB[ ,2] < lo[2]-tol  |  RGB[ ,2] > hi[2]+tol  |  RGB[ ,3] < lo[3]-tol  |  RGB[ ,3] > hi[3]+tol)

        #   clamp RGB to domain
        lo      = domain[1, ]       #lo  = -1.e-7
        hi      = domain[2, ]       #hi  = 1 + 1.e-7        

        if( length(lo) == 1 )
            {
            #   make this quick
            RGB[ RGB<lo ] = lo
            RGB[ RGB>hi ] = hi
            }
        else
            {
            #   make 2 big matrices, much slower
            n   = nrow(RGB)
            
            lo  = matrix( lo, n, 3, byrow=TRUE )
            hi  = matrix( hi, n, 3, byrow=TRUE )
            
            RGB = pmin( pmax(RGB,lo), hi )
            }

        RGBsignal  =  transfer( TF, RGB )
        }


    RGBsignal = maxSignal * RGBsignal

    colnames(RGBsignal) = c('R','G','B')

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
    
    # verify which
    w   = endIndex(which)
    if( is.na(w) )      return(NULL)            
    
    # verify maxSignal
    ok  = is.numeric(maxSignal)  &&  length(maxSignal)==1  &&  0<maxSignal
    if( ! ok )
        {
        log_string( ERROR, "maxSignal='%s' is not a positive number.", as.character(maxSignal) )
        return(NULL)
        }    
        
    rnames  = rownames(RGB)
    if( is.null(rnames)  ||  0<anyDuplicated(rnames)  )
        #   rnames is no good !  Use trivial names instead.
        rnames = 1:nrow(RGB)        
                         
    RGB = RGB/maxSignal
    

    if( is.null(TF) )
        {
        # verify space 
        idx = spaceIndex( space )
        if( is.na(idx) )    return(NULL)

        theSpace   = p.ListRGB[[idx]]

        if( w == 2 )
            TF  = theSpace$EOTF         # forwards to display
        else if( is.invertible(theSpace$OETF) )
            TF  = (theSpace$OETF)^-1    # backwards to scene
        else 
            {
            log_string( ERROR, "For RGB space='%s', the OETF is not invertible.", space )
            return(NULL)
            }    
        }
    else if( is.numeric(TF) &&  length(TF)==1  &&  0 < TF )
        {
        if( TF == 1 )
            TF  = identity.TF           # no clamping
        else if( w == 1 )
            TF  = power.OETF(TF)^-1     # backward to scene
        else
            TF  = power.EOTF(TF)        # forward to display
        }

    if( ! is.TransferFunction( TF ) )
        {
        log_string( ERROR, "argument TF is invalid."  )
        return(NULL)
        }        


    #  perhaps TODO ? -- get the domain of TF and set lo and hi from that domain
    lo  = 0
    hi  = 1
    OutOfGamut  = (RGB[ ,1] < lo  |  RGB[ ,1] > hi  |  RGB[ ,2] < lo  |  RGB[ ,2] > hi  |  RGB[ ,3] < lo  |  RGB[ ,3] > hi)
    
    if( is.identity( TF ) )
        {
        RGBlin = RGB        # TF is the identity, so copy without clamping        
        }
    else
        {
        #   clamp RGB to the unit cube (inside the normalized signal RGB gamut)
        RGB[ RGB<0 ]    = 0
        RGB[ RGB>1 ]    = 1
        
        RGBlin  = transfer( TF, RGB )   #; print( rownames(RGBlin) )        
        }
        
    
    colnames(RGBlin) = c('R','G','B')

    out = data.frame( row.names=rnames )
    out$RGB             = RGBlin
    out$OutOfGamut      = OutOfGamut
    
    return( out )
    }
    

endIndex <- function( which )
    {
    w   = pmatch( tolower(which), c('scene','display') )
    if( is.na(w) )
        {
        log_string( ERROR, "which='%s' is invalid.", which )
        return(NA_integer_)
        }

    return( w )
    }
    
    


