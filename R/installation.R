
p.ListRGB           = list()    


#   each RGB spaces is a list with these items:
#       space           the name of the space
#       primaries       4x2 matrix, with chromaticities of RGB and W
#       whiteXYZ        XYZ of the whitepoint
#       RGB2XYZ         3x3 matrix taking RGB to XYZ
#       XYZ2RGB         3x3 matrix taking XYZ to RGB
#       OETF            a number (the gamma) or a transfer function
#       OETFinv         a number (the same gamma) or a transfer function which is the inverse of OETF
#       OETFinv.fitted  present when OETFinv is a function; it is the gamma that gives the best fit to OETFinv
#       EOTF            a number (the gamma) or a transfer function
#       EOTFinv         a number (the same gamma) or a transfer function which is the inverse of EOTF
#       EOTF.fitted     present when EOTF is a function; it is the gamma that gives the best fit to EOTF
#       OOTF            a character string that describes gamma of the OOTF (system gamma), with 2 decimal places


installRGB  <-  function( space, primaries, white, OETF, EOTF=NULL, overwrite=FALSE )
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
        white.xy    = XYZ2xyY( white )[1:2]
        white.XYZ   = white
        }
        
    dim(white.XYZ)  = NULL
        

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
    
    
    primary     = primary[1:3,1:3]
    RGB2XYZ     = projectiveMatrix( t(primary), white.XYZ )
    if( is.null(RGB2XYZ) )
        {
        log.string( ERROR, "primaries is not full-rank. Please check for non-degenerate triangle with white point in interior." )
        return(FALSE)
        }
    
    rownames(RGB2XYZ)   = c('X','Y','Z')
    colnames(RGB2XYZ)   = c('R','G','B')

    theSpace$RGB2XYZ    = RGB2XYZ    
    theSpace$XYZ2RGB    = solve(RGB2XYZ) 
    
    
    
    #----   verify OETF    ----#        
    if( is.numeric(OETF)  &&  length(OETF)==1  &&  0<OETF )
        {
        theSpace$OETF       = OETF
        theSpace$OETFinv    = OETF        
        }
    else if( is.function(OETF)  )
        {
        theSpace$OETF       = OETF
        }
    else if( is.list(OETF)  &&  length(OETF)==2  &&  is.function(OETF[[1]])  &&  is.function(OETF[[2]]) )
        {
        theSpace$OETF       = OETF[[1]]
        theSpace$OETFinv    = OETF[[2]]
        }
    else if( is.character(OETF)  &&  length(OETF)==1 )
        {
        theList = functionPairFromString(OETF)
        if( is.null(theList) )  return(FALSE)
        
        theSpace$OETF       = theList$OETF
        theSpace$OETFinv    = theList$OETFinv
        }
    else
        {
        log.string( ERROR, "OETF must be a positive number, or a transfer function, or a pair of transfer functions (inverses), or one of the strings 'sRGB' etc." )
        return(FALSE)
        }

    if( is.function(theSpace$OETF)  )
        {
        #   validate OETF
        if( ! validTF(theSpace$OETF) )
            {
            log.string( ERROR, "OETF is not a valid Transfer Function" )
            return(FALSE)
            }
            
        if( is.null(theSpace$OETFinv) )
            {
            theSpace$OETFinv = makeInverseTF( theSpace$OETF )
            if( is.null(theSpace$OETFinv) ) return(FALSE)
            }
            
        #   validate OETFinv
        if( ! validTF(theSpace$OETFinv) )
            {
            log.string( ERROR, "OETFinv is not a valid Transfer Function" )
            return(FALSE)
            }

        #   validate the pair
        if( ! validTF_pair(theSpace$OETF,theSpace$OETFinv) )
            {
            log.string( ERROR, "OETF and OETFinv are valid Transfer Functions, but not inverses of each other." )
            return(FALSE)
            }
            
        theSpace$OETFinv.fitted    = fittedGammaL1( theSpace$OETFinv )
        }


    #----   verify EOTF    ----#        
    if( is.null(EOTF) )
        {
        #   use the 2 functions, or numbers, that have already been validated
        if( is.numeric(theSpace$OETF) )
            {
            theSpace$EOTF       = theSpace$OETF 
            theSpace$EOTFinv    = theSpace$OETF
            }
        else
            {
            theSpace$EOTF       = theSpace$OETFinv  # function
            theSpace$EOTFinv    = theSpace$OETF
            }
        }
    else if( is.numeric(EOTF)  &&  length(EOTF)==1  &&  0<EOTF )
        {
        theSpace$EOTF       = EOTF
        theSpace$EOTFinv    = EOTF        
        }
    else if( is.function(EOTF) )
        {
        theSpace$EOTF       = EOTF
        }
    else if( is.list(EOTF)  &&  length(EOTF)==2  &&  is.function(EOTF[[1]])  &&  is.function(EOTF[[2]]) )
        {
        theSpace$EOTF       = EOTF[[1]]
        theSpace$EOTFinv    = EOTF[[2]]
        }
    else if( is.character(EOTF)  &&  length(EOTF)==1 )
        {
        theList = functionPairFromString(EOTF)
        if( is.null(theList) )  return(FALSE)
        
        theSpace$EOTF       = theList$OETFinv
        theSpace$EOTFinv    = theList$OETF
        }
    else
        {
        log.string( ERROR, "EOTF must be a positive number, or a transfer function, or a pair of transfer functions (inverses), or one of the strings 'sRGB' etc." )
        return(FALSE)
        }
    
    if( is.function(theSpace$EOTF) )
        {
        #   validate EOTF
        if( ! validTF(theSpace$EOTF) )
            {
            log.string( ERROR, "EOTF is not a valid Transfer Function" )
            return(FALSE)
            }
            
        if( is.null(theSpace$EOTFinv) )
            {
            theSpace$EOTFinv = makeInverseTF( theSpace$EOTF )
            if( is.null(theSpace$EOTFinv) ) return(FALSE)
            }
            
        #   validate EOTFinv
        if( ! validTF(theSpace$EOTFinv) )
            {
            log.string( ERROR, "EOTFinv is not a valid Transfer Function" )
            return(FALSE)
            }

        #   validate the pair
        if( ! validTF_pair(theSpace$EOTF,theSpace$EOTFinv) )
            {
            log.string( ERROR, "EOTF and EOTFinv area Transfer Functions, but not inverses of each other." )
            return(FALSE)
            }
            
        theSpace$EOTF.fitted    = fittedGammaL1( theSpace$EOTF )            
        }
        
        
    #   OOTF is always a character string, which is simply printed in summaryRGB()
    if( is.null(EOTF) )
        ootf    = '1'  # 'identity'
    else
        {
        g1  = ifelse( is.numeric(theSpace$OETFinv), theSpace$OETFinv, theSpace$OETFinv.fitted )

        g2  = ifelse( is.numeric(theSpace$EOTF), theSpace$EOTF, theSpace$EOTF.fitted )
        
        ootf = sprintf( "%.2f", g2 / g1 )
        
        if( ! is.numeric(theSpace$OETFinv)  ||  ! is.numeric(theSpace$EOTF) )
            #  prepend a '~' which means the best-fit gamma
            ootf    = paste( '~', ootf, sep='' )
        }
    theSpace$OOTF   = ootf
        
        
    #   finally OK to install the RGB space
    p.ListRGB[[ space ]] <<- theSpace
    
    return( invisible(TRUE) )
    }
    
    
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
        log.string( ERROR, "space='%s' is unknown.", space )
        return(NULL)
        }            

    return(out)
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
    
    
getRGB  <-  function( space, full=TRUE )
    {   
    idx = spaceIndex(space)
    if( is.na(idx) )    return(NULL)

    if( ! full )
        {
        #   just the first 5 items, and no transfer functions
        return( p.ListRGB[[idx]][ 1:5 ] )
        }


    out    = p.ListRGB[[idx]]
    
    if( is.numeric(p.ListRGB[[idx]]$OETF) )
        out$OETF    = function( x ) { x^(1/p.ListRGB[[idx]]$OETF) } 
    else
        out$OETF    = p.ListRGB[[idx]]$OETF

    if( is.numeric(p.ListRGB[[idx]]$EOTF) )
        out$EOTF    = function( x ) { x^(p.ListRGB[[idx]]$EOTF) } 
    else
        out$EOTF    = p.ListRGB[[idx]]$EOTF
        
    # for OOTF, replace character string by a function
    if( out$OOTF == '1' )
        out$OOTF    = function(x) { x }     # identity
    else        
        out$OOTF    = function(x) { out$EOTF(out$OETF(x)) } # not identity


    #   erase ones we don't want
    out$OETFinv         = NULL
    out$EOTFinv         = NULL
    out$OETFinv.fitted  = NULL
    out$EOTF.fitted     = NULL


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
    
    #   add column for OETF
    myfun <- function( s )
        {
        if( is.numeric(s$OETFinv) )
            out = sprintf( "1/%.2f", s$OETFinv )
        else
            out = sprintf( "1/~%.2f", s$OETFinv.fitted )
        return( out )
        }
    out$OETF   = sapply( p.ListRGB, myfun )
    
    
    #   add column for EOTF
    myfun <- function( s )
        {
        if( is.numeric(s$EOTF) )
            out = sprintf( "%.2f", s$EOTF )
        else
            out = sprintf( "~%.2f", s$EOTF.fitted )
        return( out )
        }
    out$EOTF   = sapply( p.ListRGB, myfun )

    
    #   add column for OOTF
    myfun <- function( s )
        {
        return( s$OOTF )
        }
    out$OOTF   = sapply( p.ListRGB, myfun )

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
    