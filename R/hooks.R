

#   private global variables.  The initial 'p.' means private

#   This group is assigned during .onLoad()
#   p.microbenchmark    logical value, whether the package microbenchmark is loaded. 

#   This group is stored in system.rda
#   p.ListRGB           list of installed RGB spaces.  It must be unlocked.



p.microbenchmark    = FALSE

    
.onLoad <- function( libname, pkgname )
    {    
    #mess = '.onLoad()'
    #packageStartupMessage( mess )
    ##print( mess )
    
    
    p.microbenchmark    <<- requireNamespace( 'microbenchmark', quietly=TRUE )  #;  cat( "p.microbenchmark=", p.microbenchmark, '\n' )

    makeAllPrimaries()
    
    p.AP0_2_XYZ_MAT <<- calculateDataXYZ( AP0_PRI, 1.0)$RGB2XYZ    
    p.XYZ_2_AP0_MAT <<- calculateDataXYZ( AP0_PRI, 1.0)$XYZ2RGB

    p.AP1_2_XYZ_MAT <<- calculateDataXYZ( AP1_PRI, 1.0)$RGB2XYZ
    p.XYZ_2_AP1_MAT <<- calculateDataXYZ( AP1_PRI, 1.0)$XYZ2RGB   
    
    p.AP0_2_AP1_MAT <<- p.XYZ_2_AP1_MAT %*% p.AP0_2_XYZ_MAT
    p.AP1_2_AP0_MAT <<- p.XYZ_2_AP0_MAT %*% p.AP1_2_XYZ_MAT 
    
    p.AP1_RGB2Y     <<- p.AP1_2_XYZ_MAT[2, ]

    makeGangOf5()    
    
    makeRRTplus()
    }
    
    
.onAttach <- function( libname, pkgname )
    {
    #unlockBinding( "p.microbenchmark", asNamespace('spacesRGB') )     # asNamespace(pkgname) here generates a NOTE ! 
    #unlockBinding( "p.ListRGB", asNamespace('spacesRGB') )            # asNamespace(pkgname) here generates a NOTE !     

    #p.microbenchmark    <<- requireNamespace( 'microbenchmark', quietly=TRUE )  #;  cat( "p.microbenchmark=", p.microbenchmark, '\n' )

    #if( bindingIsLocked("p.ListRGB",asNamespace('spacesRGB')) )
    #    {
    #    packageStartupMessage( "ERROR.  Cannot unlock the list of RGB spaces.  No RGB spaces can be installed." )
    #    }
    }

    
    