
#   p.ListRGB           = list()    


#   each RGB spaces is a list with these items:
#       space           the name of the space
#
#       scene           a list with these 2 items
#           primaries   3x2 matrix, with chromaticities of RGB 
#           white       Y, xy, or XYZ of the whitepoint
#
#       display         a list with these 2 items
#           primaries   3x2 matrix, with chromaticities of RGB 
#           white       Y, xy, or XYZ of the whitepoint
#
#       OETF            a TransferFunction, with added metadata giving gamma.  negative means best fit to OETF. NA_real_ means dimension >= 2.
#       EOTF            a TransferFunction, with added metadata giving gamma.  negative means best fit to EOTF. NA_real_ means dimension >= 2.
#       OOTF            a TransferFunction, with added metadata giving gamma.  negative means best fit to OOTF. NA_real_ means dimension >= 2.


installRGB  <-  function( space, scene, display=NULL, OETF=NULL, EOTF=NULL, OOTF=NULL, overwrite=FALSE )
    {
    unlockBinding( "p.ListRGB", asNamespace('spacesRGB') )  
    
    #---   verify space ----#
    valid   = is.character(space)  &&  length(space)==1
    if( ! valid )
        {
        log_string( ERROR, "space is not a character vector of length 1." )
        return(FALSE)
        }
    
    idx = match( toupper(space), toupper(names(p.ListRGB)) )
    
    if( is.finite(idx)  &&  ! overwrite )
        {
        log_string( ERROR, "RGB space '%s' is already taken (at position %d), and overwrite is FALSE.", space, idx )
        return(FALSE)
        }
        
    ok  = is.null(scene)  ||  ( is.list(scene)  &&  length(scene)==2 )  || ( is.numeric(scene)  &&  length(scene)==8 )
    if( ! ok )
        {
        log_string( ERROR, "Argument 'scene' is not a list of length 2, or a numeric of length 8." )
        return(FALSE)
        }        
        
    ok  = is.null(display)  ||   ( is.list(display)  &&  length(display)==2 )  ||   ( is.numeric(display)  &&  length(display)==8 )
    if( ! ok )
        {
        log_string( ERROR, "Argument 'display' is not a list of length 2, or a numeric of length 8." )
        return(FALSE)
        }

    if( is.null(scene)  &&  is.null(display) )
        {
        log_string( ERROR, "Arguments 'scene' and 'display' cannot both be NULL." )
        return(FALSE)
        }
        

    theSpace   = list()    

    theSpace$space  = space   
    
    if( ! is.null(scene) )
        {
        if( is.numeric(scene) )
            scene = list( scene, 1 )
            
        theSpace$scene      = calculateDataXYZ( scene[[1]], scene[[2]] )
        
        if( is.null(theSpace$scene) )
            return(FALSE)
        }    
        
    if( ! is.null(display) )
        {
        if( is.numeric(display) )
            display = list( display, 1 )
            
        theSpace$display    = calculateDataXYZ( display[[1]], display[[2]] )
        
        if( is.null(theSpace$display) )
            return(FALSE)
        }


    
    #   go through the TF combinations
    
    if( ! is.null(OETF)  &&  ! is.null(EOTF)  &&  ! is.null(OOTF) )
        {
        log_string( ERROR, "All 3 transfer functions cannot be specified." )
        return(FALSE)
        }

    
    if( is.numeric(OETF) )
        {
        OETF    = power.OETF( OETF[1] )
        if( is.null(OETF) ) return(FALSE)
        }

    if( is.numeric(EOTF) )
        {
        EOTF    = power.EOTF( EOTF[1] )
        if( is.null(EOTF) ) return(FALSE)
        }
    else if( is.character(EOTF) )
        {
        EOTF = EOTFfromString( EOTF[1] ) 
        if( is.null(EOTF) ) return(FALSE)
        }        

        
    if( is.numeric(OOTF) )
        {
        OOTF    = power.OOTF( OOTF[1] )
        if( is.null(OOTF) ) return(FALSE)
        }
        
    listTF  = list( OETF=OETF, EOTF=EOTF, OOTF=OOTF )
    for( k in 1:length(listTF) )
        {
        tf  = listTF[[k]]
        
        if( is.null(tf) )  next
        
        if( ! is.TransferFunction(tf)  )
            {
            log_string( ERROR, "'%s' is not a valid TransferFunction.", names(listTF)[k] )
            return(FALSE)
            }
        
        if( ! (dimension(tf) %in% c(1,3))  )
            {
            log_string( ERROR, "TransferFunction '%s' is invalid, because its dimension is %d.", 
                            names(listTF)[k], dimension(tf) )
            return(FALSE)
            }
        }
        
    if( ! is.null(OETF)  &&  ! is.null(EOTF) )
        {
        OOTF    = OETF  *  EOTF
        if( is.null(OOTF) ) return(FALSE)        
        
        if( isTRUE(dimension(OOTF) == 1) )  metadata(OOTF) = NULL   # erase gamma if present
        }
    else if( ! is.null(OETF)  &&  ! is.null(OOTF) )
        {
        EOTF    = OETF^-1  *  OOTF
        if( is.null(EOTF) ) return(FALSE)        
        
        if( isTRUE(dimension(EOTF) == 1) )  metadata(EOTF) = NULL   # erase gamma if present        
        }
    else if( ! is.null(EOTF)  &&  ! is.null(OOTF) )
        {
        OETF    = OOTF  *  EOTF^-1
        if( is.null(OETF) ) return(FALSE)        
        
        if( isTRUE(dimension(OETF) == 1) )  metadata(OETF) = NULL   # erase gamma if present    
        }
    else if( ! is.null(OETF) )
        {
        EOTF    = OETF^-1
        if( is.null(EOTF) ) return(FALSE)
        metadata(EOTF) = metadata(OETF)     # for gamma
        OOTF    = identity.TF
        metadata(OOTF)  = list(gamma=1)        
        }
    else if( ! is.null(EOTF) )
        {
        OETF    = EOTF^-1
        if( is.null(OETF) ) return(FALSE)        
        metadata(OETF,add=TRUE) = metadata(EOTF)     # for gamma
        OOTF    = identity.TF
        metadata(OOTF)  = list(gamma=1)
        }
    else if( ! is.null(OOTF) )
        {
        log_string( ERROR, "The OOTF cannot be specified alone." )
        return(FALSE)
        }
    else
        {
        #   no TransferFunctions were given; set all to identity.TF
        OETF    = identity.TF
        EOTF    = identity.TF
        OOTF    = identity.TF
        
        metadata(OETF,add=TRUE)  = list(gamma=1)        
        metadata(EOTF,add=TRUE)  = list(gamma=1)        
        metadata(OOTF,add=TRUE)  = list(gamma=1)        
        }
        
    #   at this point all Transfer Functions are defined
    if( isTRUE(dimension(OETF) == 1)  &&  is.null( metadata(OETF)$gamma ) )
        {
        metadata(OETF,add=TRUE)  = list(gamma = -1 / gammaBestFit(OETF) )   # negative means best-fit
        }
        
    if( isTRUE(dimension(EOTF) == 1)  &&  is.null( metadata(EOTF)$gamma ) )
        {
        metadata(EOTF)  = list(gamma = -gammaBestFit(EOTF) )    # negative means best-fit
        }
        
    if( isTRUE(dimension(OOTF) == 1)  &&  is.null( metadata(OOTF)$gamma ) )
        {
        metadata(OOTF) = list(gamma =  -gammaBestFit(OOTF) )    # negative means best-fit
        }

        
    theSpace$OETF       = OETF
    
    theSpace$EOTF       = EOTF

    theSpace$OOTF       = OOTF

    if( is.null(display) )  
        {
        #   try to get display data from transfer function metadata
        
        for( k in 1:length(theSpace) )
            {
            tf  = theSpace[[k]]
            
            if( ! is.TransferFunction(tf) ) next
            
            primaries   = metadata(tf)$primaries      #; print(primaries)
            white       = metadata(tf)$white          #; print(white)

            if( is.null(white) )   next     # ignore this one

            if( is.null(primaries) )
                {
                #   this can happen for DCDM, the primaries are actually XYZ, and not RGB
                if( length(white) == 1 )    white = rep( white, 3 )
                
                if( length(white) != 3 )
                    {
                    log_string( ERROR, "metadata white = %s is invalid.", nicevector(white) )
                    return(FALSE)
                    }
                    
                primaries   = matrix( c(1,0,  0,1,  0,0), 3, 2, byrow=TRUE )    #;  print(primaries)
                }
            
            theSpace$display    = calculateDataXYZ( primaries, white )    #; print(  theSpace$display  )
        
            if( is.null(theSpace$display) )
                return(FALSE)        
                
            #   theSpace$display is not NULL, so we can stop searching the metadata
            break
            }
        }
    
    
    if( is.null(theSpace$display) )  theSpace$display    = theSpace$scene

    if( is.null(theSpace$scene) )    theSpace$scene      = theSpace$display

    if( is.null(theSpace$scene) )
        {
        log_string( ERROR, "Cannot assign scene primaries and white." )
        return(FALSE)        
        }
    
    if( is.null(theSpace$display) )
        {
        log_string( ERROR, "Cannot assign display primaries and white." )
        return(FALSE)        
        }
        
    #   finally OK to install the RGB space
    p.ListRGB[[ space ]] <<- theSpace
    
    return( invisible(TRUE) )
    }
    
    
    
uninstallRGB  <-  function( space )
    {   
    unlockBinding( "p.ListRGB", asNamespace('spacesRGB') )  
        
    #---   verify space ----#
    valid   = is.character(space)  &&  length(space)==1
    if( ! valid )
        {
        log_string( ERROR, "space is not a character vector of length 1." )
        return(FALSE)
        }
    
    idx = match( toupper(space), toupper(names(p.ListRGB)) )
    
    if( is.na(idx)  )
        {
        log_string( ERROR, "RGB space '%s' does not exist.", space )
        return(FALSE)
        }    
    
    #   OK to remove
    p.ListRGB[[ space ]] <<- NULL
    
    return(TRUE)
    }
    
    
getRGB  <-  function( space )       #, full=TRUE )
    {   
    idx = spaceIndex(space)
    if( is.na(idx) )    return(NULL)

    #if( ! full )
    #    #   just the first 3 items, and no transfer functions
    #    return( p.ListRGB[[idx]][ 1:3 ] )

    return( p.ListRGB[[idx]] )

    #   erase ones we don't want to return
    #out$OETF.gamma  = NULL
    #out$EOTF.gamma  = NULL
    #out$OOTF.gamma  = NULL

    #   return( out )
    }

    
getWhiteXYZ <- function( space, which='scene' )
    {
    # verify space
    idx = spaceIndex(space)
    if( is.na(idx) )    return(NULL)
    
    # verify which
    w   = endIndex(which)
    if( is.na(w) )      return(NULL)

    if( w == 1 )
        out = p.ListRGB[[idx]]$scene$whiteXYZ
    else
        out = p.ListRGB[[idx]]$display$whiteXYZ

    return( out )
    }
    
        

summaryRGB  <-  function( verbosity=1 )
    {
    if( length(p.ListRGB) == 0 )
        log_string( WARN, "There are no installed RGB spaces !" )        
        
    if( verbosity <= 0 )    return( names(p.ListRGB) )
    
    #   if( 2 <= verbosity )  print( p.ListRGB )
    
    out = data.frame( row.names=names(p.ListRGB) )

    if( length(p.ListRGB) == 0 )    return(out)
    
    
    #   add a column of chromaticities for each scene primary, including white
    pname   = rownames( p.ListRGB[[1]]$scene$primaries )
    for( i in 1:length(pname) )
        {
        mat = sapply( p.ListRGB, function(space) { space$scene$primaries[i, ] }  )     #; print(mat)
        out$X   = t(mat)
        colnames(out)[ ncol(out) ]  = sprintf( "sce%s", pname[i] )
        }
        
    #   add column of scene white XYZs
    mat = sapply( p.ListRGB, function(space) { space$scene$whiteXYZ } ) #; print(mat)
    rownames(mat) = c('X','Y','Z')
    out$sceneW  = t(mat)
    
    
    #   add a column of chromaticities for each display primary, including white
    pname   = rownames( p.ListRGB[[1]]$display$primaries )
    for( i in 1:length(pname) )
        {
        mat = sapply( p.ListRGB, function(space) { space$display$primaries[i, ] }  )     #; print(mat)
        out$X   = t(mat)
        colnames(out)[ ncol(out) ]  = sprintf( "dis%s", pname[i] )
        }
        
    #   add column of scene white XYZs
    mat = sapply( p.ListRGB, function(space) { space$display$whiteXYZ } ) #; print(mat)
    rownames(mat) = c('X','Y','Z')
    out$displayW  = t(mat)

    

    #   add column for 3 different TransferFunctions
    myfun <- function( s, tfname, fmtapprox )
        {
        n   = dimension( s[[tfname]] )
        
        if( is.na(n) )
            out = '1'       
        else if( n == 1 )
            {                
            gamma   = metadata( s[[tfname]] )$gamma
            
            if( is.null(gamma)  ||  is.na(gamma) )
                out = '1D'   
            else if( 0 < gamma )
                #   exact gamma
                out = ifelse( gamma==1, '1', sprintf( "%.2f", gamma ) )
            else
                #   approximate gamma
                out = sprintf( fmtapprox, -gamma )
            }
        else
            out = sprintf( '%dD', n )        
            
        return( out )
        }
        
        
    out$OETF   = sapply( p.ListRGB, myfun, tfname="OETF", fmtapprox="1/~%.2f" )
    
    out$EOTF   = sapply( p.ListRGB, myfun, tfname="EOTF", fmtapprox="~%.2f" )

    out$OOTF   = sapply( p.ListRGB, myfun, tfname="OOTF", fmtapprox="~%.2f" )

    return(out)
    }
    
    
#   space           name of the space
#   return value    index in p.ListRGB, partial matching and case-insensitive
#                   if no match, then NA_integer_
spaceIndex <- function( space )
    {
    ok  = is.character(space)  &&  length(space)==1
    if( ! ok )
        {
        log_string( ERROR, "space is not a character vector of length 1." )
        return(NA_integer_)
        }
        
    theNames    = names(p.ListRGB)
    if( is.null(theNames)  ||  length(theNames)==0 )
        {
        log_string( ERROR, "ERROR internal.  There are no installed RGB spaces." )
        return(NA_integer_)
        }

    idx = pmatch( toupper(space), toupper(theNames) )
    if( is.na(idx) )
        {
        log_string( ERROR, "space='%s' matches no installed spaces, or multiple spaces.", space )
        return(NA_integer_)
        }

    return(idx)
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
    
    
#############       deadwood below  #############################################
#
#   return list with 2 items
#       OETF
#       OETFinv  
functionPairFromString  <-  function( space )
    {
    out = list()
    
    if( space == 'sRGB' )
        {
        out$OETF       = OETF_sRGB
        out$OETFinv    = OETFinv_sRGB
        }
    else if( space == 'BT.709' )
        {
        out$OETF       = OETF_BT.709
        out$OETFinv    = OETFinv_BT.709
        }            
    else if( space == 'BT.2020' )
        {
        out$OETF       = OETF_BT.2020
        out$OETFinv    = OETFinv_BT.2020
        }                        
    else if( space == '240M' )
        {
        out$OETF       = OETF_240M
        out$OETFinv    = OETFinv_240M
        }                                    
    else if( space == 'ProPhotoRGB' )
        {
        out$OETF       = OETF_ProPhotoRGB
        out$OETFinv    = OETFinv_ProPhotoRGB
        }
    else
        {
        log_string( ERROR, "space='%s' is unknown.", space )
        return(NULL)
        }            

    return(out)
    }
        