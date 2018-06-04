

library( spacesRGB )

options( width=144 )
    
printf <- function( msg, ... )
    {    
    mess = sprintf( msg[1], ... )    # should this really be msg[1] ?
    cat( mess, '\n' )   #, file=stderr() )
    }
    
#   returns time in seconds, from an arbitrary origin
gettime <- function()
    {
    if( requireNamespace('microbenchmark') )
        return( microbenchmark::get_nanotime() * 1.e-9 )
    else
        return( as.double( base::Sys.time() ) )
    }
    
    
testXYZ <- function()
    {
    printf( "---------------------  testXYZ()  -----------------------" )
        
    # make random XYZs    
    set.seed(0)    
    count   = 100000
    
    RGB         = matrix( runif(3*count,max=255), ncol=3 )
    rownames(RGB)   = sprintf( "%04d", 1:count )
            
    data.space    = summaryRGB( 1 )
    
    if( nrow(data.space) == 0 )
        {
        printf(  "No RGB spaces are installed !" )    
        return(FALSE)
        }
        
    
    for( k in 1:nrow(data.space) )
        {
        space   = rownames(data.space)[k]
        
        time_start      = gettime()               
        
        
        
        XYZ             = XYZfromRGB( RGB, space=space, max=255 )$XYZ    
        RGB.back        = RGBfromXYZ( XYZ, space=space, max=255 )$RGB      #; print( 'RGB OK' )
        time_elapsed    = gettime() - time_start
        
        delta   = rowSums( abs(RGB - RGB.back) )  
        
        printf( "%s -> XYZ -> %s    max(delta)=%g   %d samples at %g sec/sample", 
                            space, space, max(delta), count, time_elapsed/count )
         
        gamma   = data.space$gamma[k]
        
        tol = ifelse( gamma=='function-pair', 5.e-12, 5.e-8 )    # pure gamma gives problems near 0 !
        
        failures = sum( tol < delta )   
        if( 0 < failures )
            {
            idx = which.max(delta)        
            printf(  "There were %d  %s -> XYZ -> %s failures.  Max error = %g",
                        failures, space, space, delta[idx] )

            df  = data.frame( row.names=1 )
            df$sRGB         = RGB[idx, ,drop=FALSE]
            df$XYZ          = XYZ[idx, ,drop=FALSE]
            df$RGB.back     = RGB.back[idx, ,drop=FALSE]        
            print( df )
            
            return(FALSE)
            }
            
        #   test pure black
        black   = c(0,0,0)
        if( ! identical( black, as.numeric(RGBfromXYZ( XYZfromRGB(black,space=space)$XYZ, space=space )$RGB ) ) )
            {
            printf(  "%s -> XYZ -> %s.back .  pure black not preserved.", space, space )
            return(FALSE)
            }        
        #   test rownames
        if( ! identical( rownames(RGB), rownames(RGB.back) ) )
            {
            printf( "%s  -> xyY -> %s .back .  rownames not preserved.", space, space )
            return(FALSE)
            }
        }

        
    return( TRUE )
    }
     
     
        
x = gettime()   # load microbenchmark
 
if( ! testXYZ() )       stop( "testXYZ() failed !", call.=FALSE )


printf(  "\nPassed all Conversion tests !" )
     