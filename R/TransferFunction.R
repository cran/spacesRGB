
#   TransferFunction
#
#   a list with items, always present
#       element     a list of elemTransferFunction's.  The length might be 0, meaning the identity function.
#       dimension   an integer
#       invertible  a logical


#   special TransferFunction, a universal identity element for all dimensions
identity.TF             <- list( element=list(), dimension=NA_integer_, invertible=TRUE )
class( identity.TF )    = c( "TransferFunction", class(identity.TF) )


#   constructor
TransferFunction  <-  function( fun, funinv, domain, range, id=NULL )
    {
    elem    = elemTransferFunction(  fun, funinv, domain, range )
    
    if( is.null(elem) )
        return(NULL)
        
    ok  = is.null(id) || (is.character(id)  &&  length(id)==1)  
    if( ! ok )
        {
        log_string( ERROR, "Argument id='%s' is not a valid character string, or NULL.", as.character(id)[1] )
        return(NULL)
        }
        
    if( is.null(id) )   id = sigfunction()
            
    #print( id )
    
    out = list()
    out$element         = list(elem)
    names(out$element)  = id
    out$dimension       = dimension(elem)
    out$invertible      = is.invertible(elem)

    class( out )    = c( "TransferFunction", class(out) )
            
    return( out )
    }
    
is.TransferFunction  <- function( x )
    {
    inherits( x, "TransferFunction" )
    }
    
as.TransferFunction.default <- function( ... )
    {
    log_string( WARN, "This function is designed to be called from other packages." )
    return(NULL)
    }    

is.identity.TransferFunction <- function(TF)
    {    
    return( length(TF$element) == 0 )
    }
    
dimension.TransferFunction <- function(TF)
    {
    return( TF$dimension )
    }    
    
is.invertible.TransferFunction <- function(TF)
    {    
    return( TF$invertible )
    }    
    
    
orientation.TransferFunction <- function(TF)
    {
    if( is.identity.TransferFunction(TF) )  return( 1 )
    
    n   = dimension(TF)
    
    theDomain   = TF$element[[1]]$domain
    
    x   = matrix( theDomain[1, ], n+1, n, byrow=TRUE )
    
    for( k in 1:n )
        x[k,k]    = theDomain[2,k]
        
    #   print(x)
        
    y   = transfer( TF, x )                                 #; print(y)
    if( is.null(y) )    return(NA_real_)
    
    mat = y[ 1:n, ] - matrix( y[n+1, ], n, n, byrow=TRUE )  #; print(mat)
    
    return( det(mat) )
    }    

domain.TransferFunction <- function(TF)
    {
    if( is.identity.TransferFunction(TF) )  return( NA_real_ )
    
    n   = dimension(TF)
    
    theDomain   = TF$element[[1]]$domain
        
    if( 1<n  &&  ncol(theDomain)==1 )
        #  replicate to m columns
        theDomain  = matrix( theDomain, 2, n )
    
    return( theDomain )
    }
    

inverse.TransferFunction <- function(TF)
    {   
    if( ! is.invertible(TF) )   return(NULL)
    
    if( is.identity(TF) )   return(TF)
    
    elements    = length( TF$element )
    
    out         = list()
    
    out$element = vector( elements, mode='list' )
    
    names( out$element )    = rev( invertNames( names(TF$element) ) )
    
    for( i in 1:elements )
        out$element[[i]]    = inverse( TF$element[[ elements+1 - i ]] )
    
    out$dimension   = TF$dimension
    out$invertible  = TF$invertible     # TRUE of course
    
    class( out )    = c( "TransferFunction", class(out) )


    return(out)
    }
    
'^.TransferFunction'  <-  function( TF, n )
    {
    if( n == 1 )    return(TF)
    
    if( n == 0 )    return( spacesRGB::identity.TF )
    
    if( n == -1 )   return( inverse(TF) )
    
    return(NULL)
    }
    
    
invertNames <- function( namevec )
    {
    for( i in 1:length(namevec) )
        {
        if( ! grepl( "\\^-1$", namevec[i] ) )
            # add decorations
            namevec[i]  = sprintf( "[%s]^-1", namevec[i] )
        else
            # remove decorations
            namevec[i]  = gsub( "(^\\[)|(\\]\\^-1$)", '', namevec[i] )
        }
    
    return( namevec )
    }
    
transfer.TransferFunction  <-  function( TF, x, domaincheck=TRUE )
    {
    if( is.identity(TF) )   return(x)
    
    if( ! is.numeric(x) )
        {
        log_string( ERROR, "argument x is not numeric.")
        return(NULL)
        }   
        
    if( length(domaincheck) == 1 )
        domaincheck = rep( domaincheck, 2 )
        
    m           = dimension(TF)

    elements    = length(TF$element)
    
    if( m == 1 )
        {
        #   every element of this transfer function is univariate, so any numeric x is acceptable
        out = as.numeric( x )
        
        for( i in 1:elements )
            {
            #   pass to TF$element[[i]]
            out = transfer( TF$element[[i]], out, domaincheck=domaincheck )
            }
                
        if( ! is.null(dim(x) ) )
            {
            dim(out)        = dim(x)    # if x is a matrix, we want output to be a matrix too, with the same dimensions of course
            rownames(out)   = rownames(x)
            }
        }
    else
        {
        #   one or more elements of this transfer function is multivariate, so dimensions must match
        x   = prepareNxM( x, m )
        if( is.null(x) )    return(NULL)
        
        out = x
        
        for( i in 1:elements )
            {
            #   pass to TF$element[[i]]
            out = transfer( TF$element[[i]], out, domaincheck=domaincheck )
            }
            
        rownames(out)   = rownames(x)    
        
        cnames  = colnames(TF$element[[elements]]$range)
        
        if( length(cnames) < m )    cnames  = rep(cnames,m)
            
        colnames(out)   = cnames 
        }
        
        
    return( out )
    }
    
#   '(.TransferFunction'  <-  function( TF, x ) { transfer(TF,x) }      # does not work !

    
    
composition.TransferFunction  <-  function( TF1, TF2 )
    {
    #   print( sys.status() )
    
    if( ! is.TransferFunction(TF2) )
        {
        log_string( ERROR, "'%s' and '%s' cannot be composed, because '%s' is not a TransferFunction object.",
                            deparse(substitute(TF1)), deparse(substitute(TF2)), deparse(substitute(TF2)) )
        return(NULL)
        }   
    
    if( is.identity(TF1) )  return(TF2)
    
    if( is.identity(TF2) )  return(TF1)
    
    #   both have a positive number of elements
    elements1   = length(TF1$element)
    elements2   = length(TF2$element)
    
    #   make a quick cancellation test that does not take long
    cancellations   = min(elements1,elements2)    #;         print( cancellations )  # optimistic
    for( i in 1:cancellations )
        {
        tf1 = TF1$element[[ elements1 - (i-1) ]]
        tf2 = TF2$element[[i]]
        
        #print( deparse( inverse(tf1)$fun ) )
        #print( deparse( tf2$fun ) )
        #print( str( inverse(tf1)$domain ) )
        #print( str( tf2$domain ) )
        #print( all.equal( inverse(tf1), tf2 ) )
        
        cancel  = is.invertible(tf1)  &&  is.invertible(tf2)  &&   identical( inverse(tf1), tf2, ignore.environment=TRUE )
        if( ! cancel )
            {
            cancellations = i-1
            break
            }
        }
        
    if( 0 < cancellations )
        {
        #   print( cancellations )
        
        #   remove cancellations items from the end of TF1
        TF1$element = TF1$element[ 1:(elements1-cancellations) ]
        elements1   = elements1 - cancellations
        
        if( 0 < elements1 )
            {
            #   recompute dimension and invertible
            TF1$dimension   = max( sapply( TF1$element, dimension.elemTransferFunction ) )
            TF1$invertible  = all( sapply( TF1$element, is.invertible.elemTransferFunction ) )
            }        
        
        
        #   remove cancellations items from the beginning of TF2
        TF2$element = TF2$element[ (1+cancellations):elements2 ]
        elements2   = elements2 - cancellations

        if( 0 < elements2 )
            {
            #   recompute dimension and invertible
            TF2$dimension   = max( sapply( TF2$element, dimension.elemTransferFunction ) )
            TF2$invertible  = all( sapply( TF2$element, is.invertible.elemTransferFunction ) )
            }
            
        if( elements1==0  &&  elements2==0 )    return( identity.TF )   # complete cancellation
            
        if( elements1 == 0 )    return( TF2 )   # only TF2 is left
        
        if( elements2 == 0 )    return( TF1 )   # only TF1 is left
        }


    #   compare dimensions
    n1  = dimension(TF1)
    n2  = dimension(TF2)
    ok  = n1==1  ||  n2==1  ||  n1==n2

    if( ! ok )
        {
        log_string( ERROR, "TF1 and TF2 cannot be composed, because dimension(TF1)=%d != %d=dimension(TF2).", n1, n2 )
        return(NULL)
        }   

    
    #   compare range1 and domain2
    range1  = TF1$element[[elements1]]$range
    domain2 = TF2$element[[1]]$domain

    inside   = insideMask( range1, domain2 )
    if( ! all(inside) )
        {
        name1   = gsub( ' ', '', names(TF1$element)[elements1] )
        name2   = gsub( ' ', '', names(TF2$element)[1] )
        
        log_string( ERROR, "%s and %s cannot be composed, because the range of %s is not inside the domain of %s, in %d of %d dimensions.",
                            name1, name2, name1, name2, length(inside)-sum(inside), length(inside) )
        return(NULL)
        }   
    

    #   finally can combine the lists !
    out = list()
    out$element     = c( TF1$element, TF2$element )         #; print( names(out$element) )
    out$dimension   = max( TF1$dimension, TF2$dimension )
    out$invertible  = TF1$invertible  &&  TF2$invertible
    
    
    #   and combine the metadata too

    class( out )    = c( "TransferFunction", class(out) )
             
    met = c( metadata(TF1), metadata(TF2) )
    if( ! is.null(met) )
        {
        dup = duplicated( names(met) )
        if( any(dup) )
            log_string( WARN, "The metadata has %d duplicates (e.g. '%s'); latter items are removed.", 
                                                    sum(dup), names(met)[ which(dup)[1] ] )            
        metadata(out)   = met[ ! dup ]
        }

    return( out )
    }
    

#   box1    2xN matrix, or a 2x1 matrix
#   box2    2xN matrix, or a 2x1 matrix
#
#   returns a logical vector of length N, whether interval i of box1 is inside interval i of box2
insideMask  <-  function( box1, box2 )
    {
    if( ncol(box1) < ncol(box2) )   box1 = matrix( box1, 2, ncol(box2) )
    if( ncol(box2) < ncol(box1) )   box2 = matrix( box2, 2, ncol(box1) )
    
    return( box2[1, ]<=box1[1, ]  &  box1[2, ]<=box2[2, ] )
    }
    
    

'*.TransferFunction'    <- function(TF1,TF2) { composition(TF1,TF2) }       # * does not imply that the composition is commutative

'%;%.TransferFunction'  <- function(TF1,TF2) { composition(TF1,TF2) }    
    
'%X%.TransferFunction'  <- function(TF1,TF2) { composition(TF1,TF2) }    
        
'%O%.TransferFunction'  <- function(TF2,TF1) { composition(TF1,TF2) }        #   math-style, order swapped 


validate.TransferFunction  <-  function( TF, points=1300, tol=5.e-7, domain=NULL )
    {
    if( is.identity(TF) )
        {
        out = TRUE
        attr( out, 'message' ) = sprintf( "'%s' is the identity.", deparse(substitute(TF)) )
        return( out )
        }
        
    ok  = is.numeric(points)  &&  length(points)==1  &&  8<=points
    if( ! ok )
        {
        log_string( ERROR, "points='%s' is invalid.", as.character(points)[1] )
        return(NULL)
        }    
    
    m   = dimension(TF)
    
    if( ! is.null(domain) )
        {
        if( ! is.matrix(domain) )
            {
            domain  = prepareNxM( domain, 2 )
            if( is.null(domain) )   return(NULL)
        
            domain  = t(domain) # the convention is to have 2 rows, not 2 columns
            
            if( 1<m  &&  ncol(domain)==1 )
                #  replicate to m columns
                domain  = matrix( domain, 2, m )            
            }
            
            
        if( ! all( dim(domain) == c(2,m) ) )
            {
            log_string( ERROR, "User-supplied domain is invalid." )
            return(NULL)
            }
            
        log_string( INFO, "Using user-supplied domain for validation." )
        }
    
    elements    = length( TF$element )  
    
    out     = logical( elements )
    mess    = vector( elements, mode='list' )
    
    for( i in 1:elements )
        {
        if( i == 1 )
            theDomain = domain
        else
            theDomain = NULL
        
        res = validate( TF$element[[i]], points=points, tol=tol, domain=theDomain )
        
        out[i]      = res
        mess[[i]]   = attr(res,'message')
        }

    attr( out, 'message' ) = mess
    
    return(out)
    }
    
print.TransferFunction  <-  function( x, ... )
    {
    #   print( sys.status() )

    if( is.identity(x) )
        {
        #theName = deparse( substitute(x) )[1]  ; print(theName)
        #if( theName == 'x' )    theName = "This"

        theName = "This"
        
        cat( sprintf( "%s is a universal identity TransferFunction.\n", theName ) )
        return( invisible(TRUE) )
        }
        
    elements    = length( x$element )
    
    for( i in 1:elements )
        {
        id  = gsub( ' ', '', names(x$element)[i] )
        cat( sprintf( "#--------------------   %s    ---------------------#\n", id ) )
        print( x$element[[i]], id )
        }
    
    return( invisible(TRUE) )
    }
    
plot.TransferFunction  <-  function( x, color='red', main=TRUE, add=FALSE, ... )
    {
    theName = nameof( x )
    
    TF  = x

    if( is.identity(TF) )
        {
        fig = NULL 
        if( add )   fig = par('fig')
        
        if( is.null(fig) )
            {
            mess    = sprintf( "'%s' is the universal identity TransferFunction.\n", theName )
            cat(mess)
            return( FALSE )
            }
            
        lines( fig[1:2], fig[3:4], col=color )
        return( invisible(TRUE) )
        }
        
    m   = dimension(TF)
    if( 4 <= m )
        {
        log_string( WARN, "'%s' has dimension %d >= 4, and cannot be plotted.", theName, m )
        return(FALSE)
        }
            
    elements    = length(TF$element)
    
    vararg  = list(...)
    log     = vararg$log
    if( is.null(log) )  log=''        
    
    theDomain   = TF$element[[1]]$domain
    theRange    = TF$element[[elements]]$range
    
    if( m == 1 )
        {
        if( ! add )
            {
            xlim    = range( theDomain )           #, 0 )
            ylim    = range( theRange )     #, 0 )
                
            plot( xlim, ylim, type='n', las=1, xlab=colnames(theDomain), ylab=colnames(theRange), log=log )
            grid( lty=1 )
            abline( h=0, v=0 )
            
            if( is.logical(main)  &&  main )
                main    = theName
                
            if( is.character(main) && 1<=length(main) )
                title( main=main[1], cex.main=1 )  
            }
        
        x   = seq( theDomain[1], theDomain[2], len=201 )
        y   = transfer( TF, x )
        lines( x, y, col=color )
        }
    else if( m == 2 )
        {
        if( ! add )
            {        
            xlim    = theRange[ ,1]
            ylim    = theRange[ ,2]
            
            plot( xlim, ylim, type='n', las=1, xlab=colnames(theRange)[1], ylab=colnames(theRange)[2], log=log )
            grid( lty=1 )
            #   abline( h=0, v=0 )
            
            if( is.logical(main)  &&  main )
                main    = theName
                
            if( is.character(main) && 1<=length(main) )
                title( main=main[1], cex.main=1 )  
            }
            
            
        #   draw approx. horizontal lines
        x1  = seq(theDomain[1,1],theDomain[2,1],len=101)
        for( x2 in seq(theDomain[1,2],theDomain[2,2],len=11) )
            {
            x1x2    = cbind( x1, x2 )
            
            y1y2    = transfer( TF, x1x2 )
            
            lines( y1y2[ ,1], y1y2[ ,2], col=color )
            }
        
        #   draw approx. vertical lines
        x2  = seq(theDomain[1,2],theDomain[2,2],len=101)
        for( x1 in seq(theDomain[1,1],theDomain[2,1],len=11) )
            {
            x1x2    = cbind( x1, x2 )
            
            y1y2    = transfer( TF, x1x2 )
            
            lines( y1y2[ ,1], y1y2[ ,2], col=color )
            }
        }
    else if( m == 3 )
        {
        if( ! requireNamespace( 'rgl', quietly=TRUE ) )
            {    
            log_string( ERROR, "Cannot plot %s, because '%s' cannot be loaded.", theName, 'rgl' )
            return( NULL )
            }

        pnts    = 11    # points on a side
        grid    = vector( 3, mode='list' )
        for( k in 1:m )
            grid[[k]]   = seq( theDomain[1,k], theDomain[2,k], length.out=pnts )
        
        x   = as.matrix( expand.grid( grid ) )  #; print( str(x) )
        
        y   = transfer( TF, x )                 #; print( str(y) )
        
        box = apply( y, 2, range )              #; print( box )
        
        dim(y)  = c( rep(pnts,3), 3 )           #; print( str(y) )
        
        #   start 3D drawing
        rgl::clear3d()        
        rgl::bg3d("gray50")
        rgl::light3d()
        
        edge    = box[2, ] - box[1, ]
        
        center  = 0.5 * (box[1, ] + box[2, ])
            
        cube    = rgl::scale3d( rgl::cube3d(col="black"), center[1], center[2], center[3] )
        cube    = rgl::translate3d( cube, center[1], center[2], center[3]  )
        rgl::wire3d( cube, lit=FALSE )

        cex     = 2
        label   = sprintf( "min=%s", nicevector(box[1, ]) )
        rgl::text3d( box[1,1], box[1,2], box[1,3], label, col='black', cex=cex, adj=1 )
        label   = sprintf( "max=%s", nicevector(box[2, ]) )        
        rgl::text3d( box[2,1], box[2,2], box[2,3], label, col='black', cex=cex, adj=0 )
        
        
        label   = colnames(theRange)
        if( ! is.null(label) )
            {
            cex     = 2
            rgl::text3d( center[1], box[1,2], box[1,3], label[1], col='black', cex=cex )
            rgl::text3d( box[1,1], center[2], box[1,3], label[2], col='black', cex=cex )
            rgl::text3d( box[1,1], box[1,2], center[3], label[3], col='black', cex=cex )        
            }
            
        for( i in 1:pnts )
            {
            for( j in 1:pnts )
                {            
                rgl::lines3d( y[i,j, ,1], y[i,j, ,2], y[i,j, ,3], col='white' )
                rgl::lines3d( y[i, ,j,1], y[i, ,j,2], y[i, ,j,3], col='white' )
                rgl::lines3d( y[ ,i,j,1], y[ ,i,j,2], y[ ,i,j,3], col='white' )                
                }
            }
        
        }
        
        
    return( invisible(TRUE) )
    }        
    
    
nameof.TransferFunction  <-  function( TF )
    {
    out = names( TF$element )
    
    out = gsub( ' ', '', out )
    
    out = paste0( out, collapse=" * " )
    
    return( out )
    }
    
    
    
    
gammaBestFit.TransferFunction <- function( TF )    
    {
    if( is.identity(TF) )   return(1)   
    
    if( dimension(TF) != 1  ) return(NA_real_)
        
    elements    = length(TF$element)   
    theDomain   = TF$element[[1]]$domain
    theRange    = TF$element[[elements]]$range
    
    ok  = all( theDomain == c(0,1) )  &&  all( theRange == c(0,1) )
    
    if( ! ok )  return(NA_real_)
        
    pnts = 512
    x   = seq( 0, 1, by=1/pnts )
    ytf = transfer( TF, x )
    
    myfun <- function( gamma )
        {
        return( sum( abs(x^gamma - ytf) ) )
        }    
    
    #   make initial estimate of gamma
    gamma   = log( ytf[pnts/2] ) / log(0.5)
        
    res = try( stats::optimize( myfun, lower=0.5*gamma, upper=2*gamma ),  silent=FALSE )   #; print( str(res) )
        
    if( inherits(res,"try-error") )    # class(res) == "try-error"    
        {
        cat( 'stats::optimize()  res = ', utils::str(res), '\n', file=stderr() )
        return( NA_real_ )
        }

    #   log_string( INFO, "LM polished gamma=%g to %g in %d iterations.", gamma, res$par, res$niter )

    return( res$minimum )    
    }
    
    
metadata.TransferFunction <- function( x, ... )
    {
    dots <- c(...) 
    
    metadata    = attr( x, "metadata" ) #   ; print(metadata)
    
    if( length(dots) == 0L )
        #   return the whole list
        metadata
    else if( length (dots) == 1L )
        #   return single item 
        metadata[dots][[1]]
    else
        #   return sublist
        metadata[dots]
    }    
    

#   value should be a named list
"metadata<-.TransferFunction" <- function( x, add=FALSE, value )
    {
    #   log.object( DEBUG, value)
    if( is.null(value) )
        {
        if( ! add ) attr(x,"metadata") <- NULL  # erase the metadata
        
        return(x) 
        }
        
        
    mask    = base::nzchar( names(value) )

    if( ! any(mask) )
        {
        print(value)
        log_string( ERROR, "none of the items in the list have names." )
        return(x)
        }

    if( ! all (mask) )
        {
        log_string( WARN, "options without name are discarded: %d", which(!mask) )
        value   = value[mask]
        }
        
        
    if( add )
        {
        metadata    = attr( x, "metadata" ) #; log.object( DEBUG, metadata )
        
        if( is.null(metadata) )
            attr(x,"metadata")  = value
        else
            attr(x,"metadata")  <- modifyList( metadata, value )
        }
    else
        {
        attr(x,"metadata") <- value 
        }
        
    return( x )
    }
    
        
    
    
#--------       UseMethod() calls           --------------#    
    
is.identity <- function(TF) 
    {
    UseMethod("is.identity")
    }    
    
is.invertible <- function(TF) 
    {
    UseMethod("is.invertible")
    }   

dimension <- function(TF) 
    {
    UseMethod("dimension")
    }  
    
nameof <- function(TF) 
    {
    UseMethod("nameof")
    }
    
orientation <- function(TF) 
    {
    UseMethod("orientation")
    }        
    
domain <- function(TF) 
    {
    UseMethod("domain")
    }            
    
validate  <-  function( TF, points=1300, tol=5.e-7, domain=NULL )    
    {
    UseMethod("validate")
    }    
    
inverse <- function(TF) 
    {
    UseMethod("inverse")
    }    
    
transfer <- function(TF,x,domaincheck=TRUE) 
    {
    UseMethod("transfer")
    }    
        
composition <- function(TF1,TF2) 
    {
    UseMethod("composition")
    } 
    
    
gammaBestFit <- function(TF)     
    {
    UseMethod("gammaBestFit")
    } 
    
metadata <- function(x,...)
    {
    UseMethod("metadata")
    }
    
"metadata<-" <- function(x,add=FALSE,value)
    {
    UseMethod("metadata<-")
    }

as.TransferFunction <- function(...)
    {
    UseMethod("as.TransferFunction")
    }
    
    
#   '%+%' <- function(TF1,TF2)  { UseMethod('%+%') }   


'%;%' <- function(TF1,TF2)  
    {
    UseMethod('%;%')
    }          
    
'%X%' <- function(TF1,TF2)  
    {
    UseMethod('%X%')
    }              
    
'%O%' <- function(TF2,TF1)
    {
    UseMethod('%O%')
    }       
    

###############     deadwood below      ########################

#    if( is.null(TF1$element) )
#        {
#        name1   = 'TF1'
#        for( where in 1:length(sys.frames()) )
#            {
#            test = paste0( deparse( substitute(TF1,sys.frame(where)) ), collapse='' )   # ; print( test )
#            if( test != name1 )
#                {
#                name1 = test
#                break
#                }
#            }
#        #print( name1 )
#        }
#    else
#        name1   = TF1$element
        