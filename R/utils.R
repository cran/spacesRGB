

#   returns time in seconds, from an arbitrary origin
gettime <- function()
    {
    if( p.microbenchmark )
        return( microbenchmark::get_nanotime() * 1.e-9 )
    else
        return( as.double( base::Sys.time() ) )
    }


#   primaries   a 3x2 matrix with RGB xy's, OR a 4x2 matrix with RGBW xy's
#   white       white xy or XYZ if primaries is 3x2, or just Y if primaries is 2x4
#
#   returns a list with 3 items:
#       primaries   4x2 matrix with all 4 chromaticities and names
#       whiteXYZ    3-vector with all names
#       RGB2XYZ     3x3 matrix
#       XYZ2RGB     3x3 matrix
#   returns NULL in case of error

calculateDataXYZ <- function( primaries, white )
    {
    #----   verify primaries    ----#
    primaries   = prepareNxM( primaries, 2 )
    if( is.null(primaries) )
        {
        log_level( ERROR, "primaries is not a 3x2 or 4x2 numeric matrix." )
        return(NULL)
        }

    if( nrow(primaries) == 3 )
        {
        #----   verify white    ----#
        valid   = is.numeric(white)  &&  length(white) %in% (2:3)
        if( ! valid )
            {
            log_level( ERROR, "white is not a numeric 2-vector or 3-vector." )
            return(NULL)
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

        primaries   = rbind( primaries, white.xy )  #   now 4x2
        }
    else if( nrow(primaries) == 4 )
        {
        #----   verify white    ----#
        valid   = is.numeric(white)  &&  length(white)==1  &&  0<white
        if( ! valid )
            {
            log_level( ERROR, "white is not a positive number." )
            return(NULL)
            }

        white.XYZ   = xyY2XYZ( c(primaries[4, ],white) )    #; print( white.XYZ )
        }
    else
        {
        log_level( ERROR, "primaries is not a 3x2 or 4x2 numeric matrix." )
        return(NULL)
        }

    #   primaries is now 4x2

    prim    = cbind( primaries, 1 - rowSums(primaries) )
    #valid   = all( 0 <= prim )
    #if( ! valid )
    #    {
    #    log_level( ERROR, "primaries does not contain 4 valid chromaticies." )
    #    print(prim)
    #    return(NULL)
    #    }

    dim(white.XYZ)  = NULL

    out = list()

    out$primaries = primaries
    rownames(out$primaries)   = c('R','G','B','W')
    colnames(out$primaries)   = c('x','y')

    names(white.XYZ)    = c('X','Y','Z')
    out$whiteXYZ        = white.XYZ

    prim        = prim[1:3,1:3]
    out$RGB2XYZ = projectiveMatrix( t(prim), white.XYZ )

    if( is.null(out$RGB2XYZ) )
        {
        log_level( ERROR, "The 4 chromaticies are degenerate. Please check for non-degenerate triangle with white point in interior." )
        return(NULL)
        }

    colnames(out$RGB2XYZ)   = c('R','G','B')
    rownames(out$RGB2XYZ)   = c('X','Y','Z')


    out$XYZ2RGB = solve(out$RGB2XYZ)

    return(out)
    }



# transfer function wrapper for calculateDataXYZ()

XYZfromRGB.TF   <-  function( primaries, white )
    {
    dataXYZ = calculateDataXYZ( primaries, white )
    if( is.null(dataXYZ) )  return(NULL)

    fun     <- function( RGB )  { as.numeric( tcrossprod( RGB, dataXYZ$RGB2XYZ ) ) }

    funinv  <- function( XYZ )  { as.numeric( tcrossprod( XYZ, dataXYZ$XYZ2RGB ) ) }

    rgbinterval = c(-65504, 65504)

    cnames  = sprintf( "linear.%s", c('R','G','B') )
    domain  = matrix( rgbinterval, 2, 3, dimnames=list(NULL,cnames) )

    #   make all vertices of cube, in 8x3 matrix
    mat     = as.matrix( expand.grid( rgbinterval, rgbinterval, rgbinterval ) )
    for( i in 1:nrow(mat) )
        mat[i, ] = fun( mat[i, ] )

    orange              = apply( mat, 2, range ) #; print(orange)
    colnames(orange)    = sprintf( "linear.%s", c('X','Y','Z') )

    out = spacesRGB::TransferFunction( fun, funinv, domain, orange, id=sigfunction() )

    metadata(out)   = list( primaries=dataXYZ$primaries, white=dataXYZ$whiteXYZ[2] )

    return( out )
    }


#   projectiveMatrix()
#
#   .matrix     invertible matrix, for example a 3x3 matrix with columns the tristimulus coordinates of RGB primaries
#   .unit       non-zero vector.  For example the tristimulus coordinates of white.
#
#   return      square matrix  B, so that
#               B = matrix  %*%  diag(a)  <=>   each column of B is a multiple of the corresponding column in .matrix
#               B %*% 1  =  .unit.      (1 is the vector of all 1s)
#
#   so for colors, B maps RGB to XYZ
#
#    Another way to write these properties:
#        B %*% I = matrix     up to multiples of the columns
#        B %*% 1  =  .unit
#   So I and 1 are the *standard* projective basis,
#   and .matrix and .unit are a different one

projectiveMatrix  <-  function( .matrix, .unit )
    {
    a   = try( solve( .matrix, .unit ), silent=TRUE )

    if( ! is.numeric(a) ) return(NULL)

    ran = range( abs(a) )   #; print(ran)

    if( ran[1] <= 1.e-6 * ran[2] ) return(NULL)

    return( .matrix  %*%  diag(a) )
    }





###########     argument processing     ##############
#
#   A   a non-empty numeric NxM matrix, or something that can be converted to be one
#
#   returns such a matrix, or NULL in case of error
#
#   This is intended to check user-supplied A, so there is a lot of checking.
#
prepareNxM  <-  function( A, M=3, Nmin=1 )
    {
    ok  = is.numeric(A) &&  0<length(A)  &&  (length(dim(A))<=2)

    ok  = ok  &&  ifelse( is.matrix(A), ncol(A)==M, ((length(A) %% M)==0)  )

    if( ! ok )
        {
        mess    = substr( as.character(A)[1], 1, 20 )

        Aname = deparse(substitute(A))

        #   in the next call, note the assignment to .topcall,
        #   which makes log_level() print the name of the parent function, and *not*  "prepareNxM()".
        log_level( ERROR, "Argument '%s' must be a non-empty numeric Nx%d matrix (with N>=%d). %s='%s...'",
                                    Aname, M, Nmin, Aname, mess, .topcall=sys.call(-1L) )
        return(NULL)
        }

    if( ! is.matrix(A) )
        A = matrix( A, ncol=M, byrow=TRUE )

    return( A )
    }


#   parent  generation of parent to return
#
#   returns a character string, with the function name and its arguments
sigfunction <- function( parent=0 )
    {
    where   = sys.parent( parent+1 )    # add 1 because parent is relative to the calling function, and not to me.   ; print(where)

    if( where == 0 )    return( "[?]" )

    out = tryCatch( deparse( sys.call(where) ), error=function(e) "[console]" )

    #   print( str( sys.call(where) )  )

    out = paste0( out, collapse='' )

    #   change spaces to Glenn style
    #   out = gsub( ' = ', '=', out )
    out = gsub( ' ', '', out )
    out = gsub( ',', ', ', out )
    out = gsub( "\\(", "( ", out )
    out = gsub( "\\)", " )", out )

    return( out )
    }




#   .pattern    a character vector of patterns
#   .string     a vector of strings
#
#   return value: a matrix of logicals
#   a row for each pattern and a column for each string
#   value of each entry is whether the corresponding string matches the corresponding pattern
multiPatternMatch <- function( .pattern, .string, .ignore=FALSE )
    {
    out = matrix( FALSE, length(.pattern), length(.string) )

    for( i in 1:length(.pattern) )
        out[i, ]    = grepl( .pattern[i], .string, ignore.case=.ignore )

    rownames(out)   = .pattern
    colnames(out)   = .string

    return(out)
    }
