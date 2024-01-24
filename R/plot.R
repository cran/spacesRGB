
#   obj             Nx3 matrix of RGBs, with rownames and colnames assigned
#                   or a data.frame with a matrix column RGB, optional columns LEFT,TOP,WIDTH,HEIGHT
#   space           name of the RGB space
#   which           'signal', 'scene', or 'display'
#   maxColorValue   of the RGB data
#   background      color of the background
#   labels          draw rownames as labels
#   shape           'full', 'half', 'hhex', 'vhex', etc.
#   add             add to existing plot

plotPatchesRGB  <-  function( obj, space='sRGB', which='signal', maxColorValue=1,
                                background='gray50', shape='full', add=FALSE,
                                labels=FALSE, ... )
    {
    if( is.matrix(obj) )
        {
        if( ncol(obj) != 3 )
            {
            log_string( ERROR, "'%s' has %d columns, but it must have 3.", deparse(substitute(obj)), ncol(obj) )
            return(FALSE)
            }    
            
        if( ! is.null(colnames(obj)) )
            {
            #   examine the column names
            mat1    =   multiPatternMatch( c("^R","^G","^B"), colnames(obj), .ignore=T )  #; print(mat1)
            mat2    =   multiPatternMatch( c("R$","G$","B$"), colnames(obj), .ignore=T )  #; print(mat2)
            
            diagl   = diag( rep(TRUE,3) )

            ok  = all( mat1==diagl ) || all( mat2==diagl )
            #   ok  = all( toupper(colnames(obj))  ==  c('R','G','B') )
            
            if( ! ok )
                log_string( WARN, "The column names of matrix obj are '%s', and do not look like RGB.", 
                                    paste( colnames(obj), collapse=',' ) )
            }
                    
        obj   = as.data.frame.model.matrix( obj )
        colnames(obj) = 'RGB'
        }
        
    if( ! is.data.frame(obj) )
        {
        log_string( ERROR, "obj is invalid; neither matrix nor data.frame" )
        return(FALSE)
        }

    #   look for column RGB
    if( is.null(obj$RGB)  ||  is.null(dim(obj$RGB))  ||  ncol(obj$RGB) != 3 )
        {
        log_string( ERROR, "data is invalid; there is no column RGB with 3 columns" )
        return(FALSE)
        }    

    
    n = nrow(obj)
    
    #   put RGBs into a vector colvec
    ok  = is.numeric(maxColorValue)  &&  length(maxColorValue)==1  &&  0<maxColorValue
    if( ! ok )
        {
        log_string( ERROR, "maxColorValue='%s' is invalid.", as.character(maxColorValue[1]) )
        return(FALSE)
        }

    full    = c('scene','signal','display')
    idx     = pmatch( tolower(which), full )
    if( is.na(idx) )
        {
        log_string( ERROR, "which='%s' is invalid", as.character(which) )
        return(FALSE)
        }    
    which   = full[idx]
    
    if( which == 'signal' )
        colvec  = rgb( obj$RGB, maxColorValue=maxColorValue )
    else
        {
        ret = SignalRGBfromLinearRGB( obj$RGB/maxColorValue, space=space, which=which )
        if( is.null(ret) )  return(FALSE)
        
        colvec  = rgb( ret$RGB )
        }    
    

        
    #   put location data into 4 vectors - left, right, bottom, top
    #   and also compute xlim and ylim
    if( all( c("LEFT","TOP","WIDTH","HEIGHT") %in% colnames(obj) ) )
        {
        #    y increases going down
        top     = obj$TOP
        bottom  = obj$TOP + obj$HEIGHT        
                
        ylim = range( top, bottom )
        ylim = ylim[2:1]    # vertical flip
        
        aspect  = 1
        }
    else if( all( c("LEFT","BOTTOM","WIDTH","HEIGHT") %in% colnames(obj) ) )
        {
        #   y increases going up
        bottom  = obj$BOTTOM
        top     = obj$BOTTOM + obj$HEIGHT
        
        ylim = range( bottom, top )
        #   ylim = ylim[2:1]    # NO vertical flip
        
        aspect  = 1
        }        
    else
        {
        #   make vertically stacked patches on the left, with lots of room for labels on the right
        #  add extra columns LEFT,TOP,WIDTH,HEIGHT        
        obj     = cbind( obj, LEFT=0, TOP=0:(n-1), WIDTH=1, HEIGHT=1 )
        
        top     = obj$TOP
        bottom  = obj$TOP + obj$HEIGHT        
        
        #   xlim = c( 0, 3 )
        ylim = c( n, 0 )
        aspect  = NA        # this NA triggers labels on the 'rightcenter', see below
        
        #print( obj )
        #return(FALSE)
        }
        
    left    = obj$LEFT
    right   = obj$LEFT + obj$WIDTH        
        
    xlim = range( left, right )
    if( xlim[2] == 1 )
        xlim[2] = 3

        
        
    #   figure out the shape
    rectangle   = shape %in% c('full','left','right','bottom','top','half')                             # logical
    triangle    = match( shape, c("bottomright", "bottomleft", "topright", "topleft") , nomatch=0 )     # integer       
    hexagon     = grepl( "^(h|v)hex$", shape )                                                          # logical
         
    ok  = rectangle  ||  (0 < triangle)  ||  hexagon
    if( ! ok )
        {
        log_string( ERROR, "shape='%s' unknown.", shape )
        return(FALSE)
        }
    
    #par( omi=rep(0,4) )
    #par( mai=c( 0.25, 0.25, 0.25, 0.25) )   

    if( ! add )
        {
        #   check background
        if( is.numeric(background) )
            {
            if( length(background) == 1 )   background = rep( background, 3 )
            
            if( length(background) != 3 )
                {
                log_string( ERROR, "background is invalid, because length(background)==%d is not 1 or 3.", length(background) )
                return(FALSE)
                }
                
            dim(background) = c(1,3)            
            
            if( which == 'signal' )
                background  = rgb( background, maxColorValue=maxColorValue )
            else
                {
                background = SignalRGBfromLinearRGB( background/maxColorValue, space=space, which=which )
                if( is.null(background) )  return(FALSE)
                background  = rgb( background$RGB )
                }
            }
            
        if( ! is.character( background ) )
            {
            log_string( ERROR, "background is invalid" )
            return(FALSE)
            }        

        bg.prev = par( bg=background )
        
        par( mgp=c(0, 0.5, 0) )
        
        plot.new()    

        plot.window( xlim, ylim, asp=aspect )    
        }
        
        

    if( rectangle )
        {
        #   optionally move 1 or all sides
        if( shape == 'left' )
            right   = 0.5*(left + right)
        else if( shape == 'right' )
            left    = 0.5*(left + right)
        else if( shape == 'top' )
            bottom  = 0.5*(top + bottom)
        else if( shape == 'bottom' )
            top     = 0.5*(top + bottom)
        else if( shape == 'half' )
            {
            x   = 0.5*(left + right)
            y   = 0.5*(top + bottom)
            
            left    = 0.5*(left + x)
            right   = 0.5*(right + x)
            top     = 0.5*(top + y)
            bottom  = 0.5*(bottom + y)
            }
            
        rect( left, bottom, right, top, col=colvec, border=NA )
        }
    else if( 0 < triangle )
        {
        for( i in 1:n )
            {
            #   assign all 4 vertices            
            xy  = c( left[i],top[i],  right[i],top[i],  left[i],bottom[i],  right[i],bottom[i] )
            dim(xy) = c(2,4)
            
            #   now drop one column
            xy  = xy[ ,-triangle]

            polygon( xy[1, ], xy[2, ], col=colvec[i], border=NA )
            }
        }
    else if( hexagon )
        {
        matx    = matrix( c(0.5,0.5, 0,1, 0,1, 0.5,0.5, 1,0, 1,0), 6, 2, byrow=TRUE )
        maty    = matrix( c(1,0, 0.75,0.25, 0.25,0.75, 0,1, 0.25,0.75, 0.75,0.25), 6, 2, byrow=TRUE )

        if( shape == 'hhex' )
            {
            #   swap matx and maty
            mat     = matx
            matx    = maty
            maty    = mat
            }
            
        for( i in 1:n )
            {
            x   = matx %*% c(left[i],right[i])
            y   = maty %*% c(bottom[i],top[i])

            polygon( x, y, col=colvec[i], border=NA )        
            }
        }

        
    if( is.logical(labels) && labels )
        labels = ifelse( is.na(aspect), 'rightcenter', 'center' )
        
    if( is.character(labels) )
        {
        rnames  = rownames(obj)
        if( is.null(rnames) )   rnames = as.character(1:n)
        
        #   center by default
        x   = (left + right)/2
        y   = (top + bottom)/2
        
        if( grepl('right',labels,ignore.case=TRUE) )
            x = right
        else if( grepl('left',labels,ignore.case=TRUE) )
            x = left
            
        if( grepl('top',labels,ignore.case=TRUE) )
            y = top
        else if( grepl('bottom',labels,ignore.case=TRUE) )
            y = bottom
            
        text( x, y, rnames, ... )
        }
        
            
    if( ! add )
        {            
        par( bg=bg.prev )   # restore previous background       
        }
        
        
    return( invisible(TRUE) )
    }
