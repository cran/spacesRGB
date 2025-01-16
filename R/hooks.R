

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
    
    if( requireNamespace( 'logger', quietly=FALSE ) )
        {
        #   log_formatter( formatter_mine )
        #   layout_mine and appender_mine are defined in logger.R
        log_formatter( logger::formatter_sprintf, namespace=pkgname )   # force sprintf(), even if glue is installed
        log_layout( layout_mine, namespace=pkgname )                    # put fn() between timestamp and the msg    
        log_appender( appender_mine, namespace=pkgname )                # maybe stop on ERROR or FATAL
        log_threshold( WARN, namespace=pkgname )                        # default is INFO
        }    
    

    #   this one makes 7 standard primary matrices, all 4x2
    makeAllPrimaries()
    
    #   make more matrices, these are all 3x3
    p.AP0_2_XYZ_MAT <<- calculateDataXYZ( AP0_PRI, 1.0)$RGB2XYZ    
    p.XYZ_2_AP0_MAT <<- calculateDataXYZ( AP0_PRI, 1.0)$XYZ2RGB

    p.AP1_2_XYZ_MAT <<- calculateDataXYZ( AP1_PRI, 1.0)$RGB2XYZ
    p.XYZ_2_AP1_MAT <<- calculateDataXYZ( AP1_PRI, 1.0)$XYZ2RGB   
    
    p.AP0_2_AP1_MAT <<- p.XYZ_2_AP1_MAT %*% p.AP0_2_XYZ_MAT
    p.AP1_2_AP0_MAT <<- p.XYZ_2_AP0_MAT %*% p.AP1_2_XYZ_MAT 
    
    p.AP1_RGB2Y     <<- p.AP1_2_XYZ_MAT[2, ]

    #   this one makes 5 standard EOTFs; they do not require the above matrices
    makeGangOf5() 
    
    #   this one makes more matrices, transfer functions, and one more EOTF
    makeRRTplus()
    
    #   finally ready to make initial dictionary of 8 color spaces - about 20 msec
    makeInitialDictionary()    
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

    
    