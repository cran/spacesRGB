
p.ListRGB = NULL

#   an advantage of the private data in "sysdata.rda" is that these
#   do not have to be exposed, and therefore documented
savePrivateDatasets  <- function( .path="sysdata.rda" )
    {
    savevec = character(0)
    
    #white.D65   = c( 0.95047, 1, 1.08883 )
    #white.D50   = c( 0.96422, 1, 0.82521 )
        
    makeAllPrimaries()        
    makeGangOf5()
        
    p.ListRGB <<- list()

    cat( "Installing sRGB...", '\n', file=stderr() )
    #prim    = matrix( c(0.64,0.33,  0.30,0.60,  0.15,0.06 ), 3, 2, byrow=TRUE )
    #white   = c( 0.3127, 0.3290 )
    installRGB( 'sRGB', scene=REC709_PRI, EOTF=sRGB.EOTF, overwrite=TRUE )    
    
    cat( "Installing AdobeRGB...", '\n', file=stderr() )    
    prim    = c(0.64,0.33,  0.21,0.71,  0.15,0.06  )
    white   = c( 0.3127, 0.3290 )    
    installRGB( 'AdobeRGB', scene=list(prim,white), EOTF=563/256, overwrite=TRUE )

    cat( "Installing ProPhotoRGB...", '\n', file=stderr() )    
    prim    = c(0.7347,0.2653,  0.1596,0.8404,  0.0366,0.0001  )
    white   = c( 0.3457,0.3585 )    # D50
    installRGB( 'ProPhotoRGB', scene=list(prim,white), EOTF='ProPhotoRGB', overwrite=TRUE )

    cat( "Installing AppleRGB...", '\n', file=stderr() )
    prim    = c(0.625,0.34,  0.28,0.595, 0.155,0.07)
    white   = c( 0.3127, 0.3290 ) 
    installRGB( 'AppleRGB', scene=list(prim,white), EOTF=1.8, overwrite=TRUE )

    
    # this one has the same primaries as sRGB
    cat( "Installing BT.709...", '\n', file=stderr() )
    installRGB( 'BT.709', scene=REC709_PRI, EOTF=BT.709.EOTF, overwrite=TRUE )    
    
    cat( "Installing BT.2020...", '\n', file=stderr() )
    #prim    = c(0.708,0.292,  0.170,0.797,  0.131,0.046 ) 
    #white   = c( 0.3127, 0.3290 )
    installRGB( 'BT.2020', scene=REC2020_PRI, EOTF=BT.2020.EOTF, overwrite=TRUE )    
    
    cat( "Installing 240M...", '\n', file=stderr() )
    prim    = c(0.64,0.34,  0.31,0.595,  0.155,0.07 )
    white   = c( 0.3127, 0.3290 )
    installRGB( '240M', scene=list(prim,white), EOTF='240M', overwrite=TRUE )  

    # Install an RGB space named 'HD+2.4', with encoding from BT.709 and display from BT.1886.
    # the OOTF for this space is non-trivial
    cat( "Installing HD+2.4...", '\n', file=stderr() )    
    #prim    = matrix( c(0.64,0.33,  0.30,0.60,  0.15,0.06,  0.3127,0.3290), 4, 2, byrow=TRUE )
    installRGB( "HD+2.4", scene=REC709_PRI, OETF=BT.709.EOTF^-1, EOTF=2.4, overwrite=TRUE )    

    savevec = c( savevec, "p.ListRGB" ) 

    ##  finally ready to save it
    save( list=savevec, file=.path, compress='xz' )   #     'xz'  'gzip'  FALSE
    
    mess    = sprintf( "Saved the following to '%s'.", .path )
    cat( mess, '\n', file=stderr() )
    cat( savevec, '\n', file=stderr() )
        
    return( invisible(TRUE) )
    }    
    
        