

#   returns time in seconds, from an arbitrary origin
gettime <- function()
    {
    if( p.microbenchmark )
        return( microbenchmark::get_nanotime() * 1.e-9 )
    else
        return( as.double( base::Sys.time() ) )
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
prepareNxM  <-  function( A, M=3 )
    {
    ok  = is.numeric(A)  &&  ((length(A) %% M)==0)  &&  0<length(A)
    if( ! ok )
        {
        #print( "prepareNx3" )
        #print( sys.frames() )
        mess    = substr( as.character(A)[1], 1, 10 )
        #arglist = list( ERROR, "A must be a non-empty numeric Nx3 matrix (with N>0). A='%s...'", mess )
        #do.call( log.string, arglist, envir=parent.frame(n=3) )
        #myfun   = log.string
        #environment(myfun) = parent.frame(3)
        log.string( ERROR, "Argument A must be a non-empty numeric Nx%d matrix (with N>0). A='%s...'", M, mess )
        return(NULL)
        }
    
    if( is.null(dim(A)) )
        A = matrix( A, ncol=M, byrow=TRUE )
    else if( ncol(A) != M )
        A = t(A)
        
    return( A )
    }
        