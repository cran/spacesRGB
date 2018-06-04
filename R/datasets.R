


#   an advantage of the private data in "sysdata.rda" is that these
#   do not have to be documented, and therefore exposed    
savePrivateDatasets  <- function( .path="sysdata.rda" )
    {
    savevec = character(0)
    
    white.D65   = c( 0.95047, 1, 1.08883 )
    white.D50   = c( 0.96422, 1, 0.82521 )
        
    #   install 2 RGB spaces
    prim    = matrix( c(0.64,0.33,  0.30,0.60,  0.15,0.06 ), 3, 2, byrow=T )
    white   = c( 0.3127, 0.3290 )
    installRGB( 'sRGB', prim, white.D65, 'sRGB', overwrite=TRUE )    
    
    prim    = matrix( c(0.64,0.33,  0.21,0.71,  0.15,0.06  ), 3, 2, byrow=T )
    white   = c( 0.3127, 0.3290 )    
    installRGB( 'AdobeRGB', prim, white.D65, 563/256, overwrite=TRUE )
    
    
    prim    = matrix( c(0.7347,0.2653,  0.1596,0.8404,  0.0366,0.0001  ), 3, 2, byrow=T )
    white   = c( 0.3457,0.3585 )    
    installRGB( 'ProPhotoRGB', prim, white.D50, 'ProPhotoRGB', overwrite=TRUE )

    savevec = c( savevec, "p.ListRGB" ) 

    ##  finally ready to save it
    save( list=savevec, file=.path, compress='xz' )   #     'xz'  'gzip'  FALSE
    
    mess    = sprintf( "Saved the following to '%s'.", .path )
    cat( mess, '\n', file=stderr() )
    cat( savevec, '\n', file=stderr() )
        
    return( invisible(TRUE) )
    }    
    
        