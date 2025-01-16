
g.options   <- list( stoponerror = TRUE                                  #   must be logical
                    )
    
#   put fn() between timestamp and the msg    
layout_mine <- structure(
    function(level, msg, namespace="spacesRGB",
                                    .logcall = sys.call(), .topcall = sys.call(-1), .topenv = parent.frame())
        {
        # cat( "obj_addr()=", obj_addr( .topcall[[1L]] ), '\n' )
        # cat( "deparse1 =", deparse1( .topcall[[1L]] ), '\n' )
        
        fn  = deparse1( .topcall[[1L]] )
        
        paste0( attr(level, 'level'), ' [', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '] ', namespace, "::", fn, '(). ', msg )
        },
    generator = quote(layout_mine())
)



#   maybe stop on ERROR or FATAL
appender_mine <- structure(
    function(lines)
        {
        cat(lines, file = stderr(), sep = '\n' )
        
        #   test for STOP
        if( any( grepl("^(ERR|FATAL)",lines ) )  )
            {
            stop( "Stopping in package 'spacesRGB', because level is ERROR or FATAL.", call.=FALSE )
            }
        },
    generator = quote(appender_mine())
    )

