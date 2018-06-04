

#   private global variables.  The initial 'p.' means private

#   This group is programmatically created during .onAttach()
#   p.microbenchmark    logical value, whether the package microbenchmark is loaded.  It must be unlocked.
#   p.ListRGB           list of installed RGB spaces. 2 spaces added in .onAttach().  It must be unlocked.



p.microbenchmark    = FALSE

    
.onLoad <- function( libname, pkgname )
    {    
    #   unlockBinding() fails here in .onLoad(), so use .onAttach() instead

    }
    
    
.onAttach <- function( libname, pkgname )
    {
    unlockBinding( "p.microbenchmark", asNamespace('spacesRGB') )     # asNamespace(pkgname) here generates a NOTE ! 
    unlockBinding( "p.ListRGB", asNamespace('spacesRGB') )            # asNamespace(pkgname) here generates a NOTE !     

    p.microbenchmark    <<- requireNamespace( 'microbenchmark', quietly=TRUE )  #;  cat( "p.microbenchmark=", p.microbenchmark, '\n' )

    if( bindingIsLocked("p.ListRGB",asNamespace('spacesRGB')) )
        {
        packageStartupMessage( "ERROR.  Cannot unlock the list of RGB spaces.  No RGB spaces can be installed." )
        }
    
    }

    
    