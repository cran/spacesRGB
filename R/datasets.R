
# p.ListRGB = NULL

#   make p.ListRGB and save it to "sysdata.rda" 

savePrivateDatasets  <- function( .path="sysdata.rda" )
    {
    savevec = character(0)
    
    makeInitialDictionary()
    
    savevec = c( savevec, "p.ListRGB" ) 

    ##  ready to save it
    save( list=savevec, file=.path, compress='xz' )   #     'xz'  'gzip'  FALSE
    
    mess    = sprintf( "Saved the following to '%s'.", .path )
    cat( mess, '\n', file=stderr() )
    cat( savevec, '\n', file=stderr() )
        
    return( invisible(TRUE) )
    }    
    
        