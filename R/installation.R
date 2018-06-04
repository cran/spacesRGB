
p.ListRGB           = list()    


installRGB  <-  function( space, primaries, white, gamma, peakRGB=1, overwrite=FALSE )
    {
    #---   verify space ----#
    valid   = is.character(space)  &&  length(space)==1
    if( ! valid )
        {
        log.string( ERROR, "space is not a character vector of length 1." )
        return(FALSE)
        }
    
    idx = match( toupper(space), toupper(names(p.ListRGB)) )
    
    if( is.finite(idx)  &&  ! overwrite )
        {
        log.string( ERROR, "RGB space '%s' is already taken (at position %d), and overwrite is FALSE.", space, idx )
        return(FALSE)
        }
        
    theSpace   = list()        
    
    theSpace$space  = space
    
    #----   verify primaries    ----#
    valid   = is.numeric(primaries)  &&  all( dim(primaries)==c(3,2) )
    if( ! valid )
        {
        log.string( ERROR, "primaries is not a 3x2 numeric matrix." )
        return(FALSE)
        }
        
    #----   verify white    ----#
    valid   = is.numeric(white)  &&  length(white) %in% (2:3)
    if( ! valid )
        {
        log.string( ERROR, "white is not a numeric 2-vector or 3-vector." )
        return(FALSE)
        }
        
        
    if( length(white) == 2 )
        {
        white.xy    = white
        white.XYZ   = xyY2XYZ( c(white,1) )
        }
    else
        {
        if( white[2] != 1 )
            {
            log.string( ERROR, "The white Y = %g != 1.", white[2] )
            return(FALSE)
            }
        white.xy    = XYZ2xyY( white )[1:2]
        white.XYZ   = white
        }
        
    dim(white.XYZ)  = NULL
        
    #----   verify peakRGB    ----#
    valid   = is.numeric(peakRGB)  &&  length(peakRGB)==1  &&  all( 0<peakRGB  &  peakRGB<=1 )
    if( ! valid )
        {
        log.string( ERROR, "peakRGB='%s' is not a number, or is not in the interval (0,1].", as.character(peakRGB[1]) )
        return(FALSE)
        }
        
    peakRGB = rep( peakRGB, 3 )
    names(peakRGB)  = c('R','G','B')

    primaries   = rbind( primaries, white.xy )
    
    primary = cbind( primaries, 1 - rowSums(primaries) )
    valid   = all( 0 <= primary )
    if( ! valid )
        {
        log.string( ERROR, "primaries does not contain 4 valid chromaticies." )
        return(FALSE)
        }
        
    theSpace$primaries = primaries     # only the white is really required for later use
    rownames(theSpace$primaries)   = c('R','G','B','W')
    colnames(theSpace$primaries)   = c('x','y')
          
    names(white.XYZ)    = c('X','Y','Z')
    theSpace$whiteXYZ   = white.XYZ    
    
    theSpace$peakRGB    = peakRGB
    
    primary     = primary[1:3,1:3]
    RGB2XYZ     = projectiveMatrix( t(primary), white.XYZ/peakRGB )
    if( is.null(RGB2XYZ) )
        {
        log.string( ERROR, "primaries is not full-rank. Please check for non-degenerate triangle with white point in interior." )
        return(FALSE)
        }
    
    rownames(RGB2XYZ)   = c('X','Y','Z')
    colnames(RGB2XYZ)   = c('R','G','B')

    theSpace$RGB2XYZ    = RGB2XYZ    
    theSpace$XYZ2RGB    = solve(RGB2XYZ) 
    
    #----   verify gamma    ----#    
    if( is.numeric(gamma)  &&  length(gamma)==1  &&  0<gamma )
        {
        theSpace$gamma = gamma
        }
    else if( is.list(gamma) &&  length(gamma)==2  &&  is.function(gamma[[1]])  &&  is.function(gamma[[2]]) )
        {
        #   check that 0->0 and 1->1
        endpoint    = c(0,1)
        for( k in 1:2 )
            {
            ok  = all( gamma[[k]](endpoint) == endpoint )
            if( ! ok )
                {
                log.string( ERROR, "gamma transfer function %d does not map 0->0 and 1->1.", k )
                return(FALSE)
                }
                
            #   check that functions preserve matrix dimensions                
            ok  = all( dim(gamma[[k]](primaries)) == dim(primaries) )
            if( ! ok )
                {
                log.string( ERROR, "gamma transfer function %d does not preserve matrix dimensions.", k )
                return(FALSE)
                }
            }

        theSpace$gamma = gamma
        }
    else if( is.character(gamma)  &&  length(gamma)==1 )
        {
        if( gamma == 'sRGB' )
            #   install the special sRGB functions.  No need to check them.
            theSpace$gamma = list( LinearFromDisplay_sRGB, DisplayFromLinear_sRGB )
        else if( gamma == 'ProPhotoRGB' )
            #   install the special ProPhotoRGB functions.  No need to check them.
            theSpace$gamma = list( LinearFromDisplay_ProPhotoRGB, DisplayFromLinear_ProPhotoRGB )
        else
            {
            log.string( ERROR, "gamma='%s' is invalid.", gamma )
            return(FALSE)
            }            
        }
    else
        {
        log.string( ERROR, "gamma must be a positive number, or a pair of valid transfer functions, or one of the strings 'sRGB' or 'ProPhotoRGB'." )
        return(FALSE)
        }

    #   finally OK to install the RGB space
    p.ListRGB[[ space ]] <<- theSpace
    
    return(TRUE)
    }
    
    
uninstallRGB  <-  function( space )
    {   
    #---   verify space ----#
    valid   = is.character(space)  &&  length(space)==1
    if( ! valid )
        {
        log.string( ERROR, "space is not a character vector of length 1." )
        return(FALSE)
        }
    
    idx = match( toupper(space), toupper(names(p.ListRGB)) )
    
    if( is.na(idx)  )
        {
        log.string( ERROR, "RGB space '%s' does not exist.", space )
        return(FALSE)
        }    
    
    #   OK to remove
    p.ListRGB[[ space ]] <<- NULL
    
    return(TRUE)
    }
    
    
getRGB  <-  function( space )
    {   
    idx = spaceIndex(space)
    if( is.na(idx) )    return(NULL)

    #if( short )
    #    return( list( space=p.ListRGB[[idx]]$space, primaries=p.ListRGB[[idx]]$primaries  ) )

    out    = p.ListRGB[[idx]]
    
    if( is.numeric(p.ListRGB[[idx]]$gamma) )
        {
        out$EOCF    = function( x ) { x^(p.ListRGB[[idx]]$gamma) }        
        out$OECF    = function( x ) { x^(1/p.ListRGB[[idx]]$gamma) }
        }
    else
        {
        out$EOCF    = p.ListRGB[[idx]]$gamma[[1]]        
        out$OECF    = p.ListRGB[[idx]]$gamma[[2]]
        }

    out$gamma = NULL

    return( out )
    }

    
    
summaryRGB  <-  function( verbosity=1 )
    {
    if( length(p.ListRGB) == 0 )
        log.string( WARN, "There are no installed RGB spaces !" )        
        
    if( verbosity <= 0 )    return( names(p.ListRGB) )
    
    #   if( 2 <= verbosity )  print( p.ListRGB )
    
    out = data.frame( row.names=names(p.ListRGB) )

    if( length(p.ListRGB) == 0 )    return(out)
    
    
    #   add a column of chromaticities for each primary, including white
    pname   = rownames( p.ListRGB[[1]]$primaries )
    
    for( i in 1:length(pname) )
        {
        mat = sapply( p.ListRGB, function(s) { s$primaries[i, ] }  )     #; print(mat)
    
        out$X   = t(mat)
        
        colnames(out)[ ncol(out) ]  = pname[i]
        }
        
    #   add column of white XYZs
    mat = sapply( p.ListRGB, function(s) { s$whiteXYZ } )
    rownames(mat) = c('X','Y','Z')
    out$white   = t(mat)
    
    out$gamma   = sapply( p.ListRGB, function(s) { ifelse( is.numeric(s$gamma),as.character(s$gamma),'function-pair') } )
    
    return(out)
    }

    
#############     sRGB  EOCF and OECF [0,1]  <-->  [0,1]   ###############    
#   maps [0,1] to [0,1].   Also extended to negative numbers in a symmetric way
#   it is OK if input is a matrix, and then the return value is a matrix of the same shape
DisplayFromLinear_sRGB <- function( lin )
    {
    s = sign(lin)    
    out = s * lin    #   now non-negative
    out = ifelse( out <= 0.0031308,    12.92 * out,  (1055 * out^(1/2.4) - 55 ) / 1000 )
    return( s * out )
    }

LinearFromDisplay_sRGB <- function( disp )
    {
    s = sign(disp)    
    out = s * disp    #   now non-negative
    out = ifelse( out <= 12.92*0.0031308,   out / 12.92,   ( (1000*out + 55)/(1055) ) ^ 2.4  )
    return( s * out )
    }    
    


#############     ProPhotoRGB  EOCF and OECF [0,1]  <-->  [0,1]   ###############    
#   maps [0,1] to [0,1].   Also extended to negative numbers in a symmetric way
#   it is OK if input is a matrix, and then the return value is a matrix of the same shape
DisplayFromLinear_ProPhotoRGB <- function( lin )
    {
    s = sign(lin)    
    out = s * lin    #   now non-negative
    out = ifelse( out <= 1/512,  16 * out,  out^(1/1.8) )
    return( s * out )
    }

LinearFromDisplay_ProPhotoRGB <- function( disp )
    {
    s = sign(disp)    
    out = s * disp    #   now non-negative
    out = ifelse( out <= 1/32,   out / 16,  out ^ 1.8  )
    return( s * out )
    }    
    
    
    


xyY2XYZ <- function( xyY )
    {
    xyY = prepareNxM(xyY)
    if( is.null(xyY) )  return(NULL)

    XYZ <- cbind( NA_real_, xyY[,3], NA_real_)
    rownames(XYZ) = rownames(xyY)
    colnames(XYZ) = c('X','Y','Z')

    w <- which( xyY[,2] != 0  &  0 <= xyY[,3] )    # was 0 < xyY[ ,2]

    if (length(w)>0){
    xyY_sub = xyY[w, ,drop=FALSE]
    mult    =  xyY_sub[ ,3] / xyY_sub[ ,2]
    XYZ[w,1] <- mult * xyY_sub[ ,1]
    XYZ[w,3] <- mult * (1-xyY_sub[ ,1]-xyY_sub[ ,2])
    }

    #   treat Y=0 as a special case - pure black
    w <- which( xyY[,3] == 0 )
    if( length(w) > 0 )
        XYZ[w,1:3] = 0

    XYZ
    }
    
XYZ2xyY <- function( XYZ )
    {
    XYZ = prepareNxM(XYZ)
    if( is.null(XYZ) )  return(NULL)
    
    xyY <- cbind(NA_real_, NA_real_, XYZ[ ,2])
    rownames(xyY) = rownames(XYZ)
    colnames(xyY) = c('x','y','Y')

    denom       = rowSums( XYZ )
    w <- which(0<denom  &  0<=XYZ[ ,2])
    if (length(w)>0){
    xyY[w,1]    = XYZ[w,1] / denom[w]    
    xyY[w,2]    = XYZ[w,2] / denom[w]    
    }
    
    xyY
    }
    