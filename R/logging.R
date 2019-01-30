
FATAL <- 1L
names(FATAL) <- "FATAL"
ERROR <- 2L
names(ERROR) <- "ERROR"
WARN <- 4L
names(WARN) <- "WARN "
INFO <- 6L
names(INFO) <- "INFO "
DEBUG <- 8L
names(DEBUG) <- "DEBUG"
TRACE <- 9L
names(TRACE) <- "TRACE"

log.string <- function( level, msg, ... )
    {    
    conn    = stderr()
    
    if( ! is.integer(level) )
        {
        cat( "log.string(). level is not an integer.\n", file=conn )
        return( invisible(FALSE) )
        }
        
    generation  = 1L
    if( 2 <= length(level) )
        {
        #   hack to get higher generation parents !!
        generation  = level[2]
        level       = level[1]  # this preserves names(level[1])
        }          
        
    msg = sprintf( msg[1], ... )    # should this really be msg[1] ?

    #   find the name of the calling function
    #print( "log.string" )
    # print( sys.status() )
        
    where   = sys.parent(generation)  # ; print(where)
  
    if( 0 < where )
        namefun = tryCatch( deparse(sys.call(where)[[1L]]), error=function(e) "[console]" )
    else
        namefun = "[??]"

    if( ! grepl( "^spacesRGB", namefun ) )
        namefun = paste0( "spacesRGB::", namefun, collapse='' )
        
    mess    = paste0( namefun, "(). ", names(level), ".  ", msg, collapse='' )

    if( level <= FATAL )
        stop( mess, '\n',  call.=FALSE )

    cat( mess, '\n', file=conn ) ;   flush(conn)

    return( invisible(TRUE) )
    }
