
#   class elemTransferFunction
#
#   these are not designed to be user-callable; they are not exported
#
#   a list with items
#       fun     always present, never trivial
#       funinv  possibly missing
#       domain  a 2xN matrix.  all must be finite; Inf, NA, NaN not allowed
#       range   a 2xN matrix


#   constructor
elemTransferFunction  <-  function( fun, funinv, domain, range )
    {
    if( ! is.function(fun) )
        {
        log_string( ERROR, "Argument fun is not a valid function." )
        return(NULL)
        }

    ok  = is.null(funinv) || is.function(funinv)
    if( ! ok )
        {
        log_string( ERROR, "Argument fun is not a valid function, or NULL." )
        return(NULL)
        }

    domain  = prepareBox(domain)
    if( is.null(domain) )   return(NULL)


    range   = prepareBox(range)
    if( is.null(range) )   return(NULL)


    n   = ncol(domain)
    if( n != ncol(range))
        {
        log_string( ERROR, "dimension of domain = %d  !=  %d = dimension of range.", n, ncol(range) )
        return(NULL)
        }

    if( n == 1 )
        {
        #   check monotonicity
        x       = seq( domain[1], domain[2], length.out=257 )
        y       = fun(x)

        if( ! all(is.finite(y) ) )
            {
            idx = which( ! is.finite(y) )[1]
            log_string( ERROR, "x=%g  maps to y=%g, which is not finite.", x[idx], y[idx] )
            return(NULL)
            }

        delta   = diff(y)
        ok  = all(delta<0)  ||  all(0<delta)
        if( ! ok )
            {
            log_string( ERROR, "The univariate function is not monotone."  )
            return(NULL)
            }

        if( is.null(funinv) )
            {
            #   make an approximate inverse
            funinv  = makeInverseTF( fun, domain )
            log_string( INFO, "Created an approximate inverse for '%s', using stats::splinefun().", deparse(substitute(fun)) )
            }
        }


    if( is.null(colnames(domain)) )
        {
        #   supply some good defaults
        if( n == 1 )
            colnames(domain)    = 'x'
        else
            colnames(domain)    = paste( 'x', 1:n, sep='.' )
        }

    if( is.null(colnames(range)) )
        {
        #   supply some good defaults
        if( n == 1 )
            colnames(range) = 'y'
        else
            colnames(range) = paste( 'y', 1:n, sep='.' )
        }


    out = list()
    out$fun     = fun
    out$funinv  = funinv
    out$domain  = domain
    out$range   = range

    class( out )    = c( "elemTransferFunction", class(out) )

    return( out )
    }


#   box         a 2xN matrix, or a 2-vector
#
#   returns a 2xN matrix
prepareBox <- function( box )
    {
    if( ! is.numeric(box) )
        {
        #   notice hack to make log_string() print name of parent function
        log_string( c(ERROR,2L), "box '%s' is invalid, because it is not numeric.", deparse(substitute(box)) )
        return(NULL)
        }

    if( ! all( is.finite(box) ) )
        {
        log_string( ERROR, "The box is invalid, because some values are not finite."  )
        return(NULL)
        }


    if( is.matrix(box)  &&  nrow(box)==2  &&  0<ncol(box) )
        {
        # good to go, do nothing
        }
    else if( length(box) == 2 )
        {
        #   box = matrix( box, nrow=2, byrow=FALSE )
        dim(box)    = c(2,1)
        }
    else
        {
        #   notice hack to make log_string() print name of parent function
        log_string( c(ERROR,2L), "box '%s' is not valid.  It must be a 2xM matrix, or a 2-vector.", deparse(substitute(box)) )
        return(NULL)
        }

    #   check orientation of all the intervals
    len = box[2, ] - box[1, ]

    valid   = 0 < len
    if( ! all(valid) )
        {
        #   notice hack to make log_string() print name of parent function
        log_string( c(ERROR,2L), "box '%s' is not valid. %d of %d intervals have incorrect order.",
                        deparse(substitute(box)), length(valid)-sum(valid), length(valid) )
        return(NULL)
        }

    rownames(box)   = c('min','max')

    return(box)
    }


dimension.elemTransferFunction <- function(TF)
    {
    return( ncol(TF$domain) )
    }

is.invertible.elemTransferFunction <- function(TF)
    {
    return( ! is.null(TF$funinv) )
    }


inverse.elemTransferFunction <- function(TF)
    {
    if( ! is.invertible(TF) )   return(NULL)

    #out = list()
    #out$fun     = TF$funinv
    #out$funinv  = TF$fun
    #out$domain  = TF$range
    #out$range   = TF$domain

    out = list( fun=TF$funinv, funinv=TF$fun, domain=TF$range, range=TF$domain )

    class( out )    = c( "elemTransferFunction", class(out) )

    return(out)
    }


#   x   an NxM matrix.  For max speed no reshaping is done, since this is taken care of by the caller.

transfer.elemTransferFunction  <-  function( TF, x, domaincheck=c(TRUE,TRUE) )
    {
    #   if( is.identity(TF) )   return(x)

    if( ! is.numeric(x) )
        {
        log_string( ERROR, "argument x is not numeric.")
        return(NULL)
        }

    if( length(domaincheck) != 2 )
        domaincheck = rep( domaincheck[1], 2 )

    m   = dimension(TF)

    if( m == 1 )
        {
        #   this transfer function is univariate, so any numeric x is acceptable

        if( domaincheck[1] )
            #   set values below lower limit to NA
            x[ x<TF$domain[1] ] = NA

        if( domaincheck[2] )
            #   set values above upper limit to NA
            x[ TF$domain[2]<x ] = NA

        #   ready to pass to fun
        out = TF$fun(x)

        dim(out)    = dim(x)    # if x is a matrix, we want output to be a matrix too, with the same dimensions of course

        return( out )
        }

    #   domain is multivariate, so dimensions must match
    ok  = is.matrix(x)  &&  ncol(x) == m
    if( ! ok )
        {
        log_string( ERROR, "argument x is not an Nx%d matrix.", m )
        return(NULL)
        }


    #   transform one row of x at a time
    out = matrix( NA_real_, nrow(x), ncol(x) )
    rownames(out)   = rownames(x)
    colnames(out)   = colnames(TF$range)

    for( i in 1:nrow(x) )
        {
        xi  = x[i, ]    # extract row i

        if( any( is.na(xi) ) )  next

        #   optionally check that xi is in the domain box
        #   print( domaincheck )

        if( domaincheck[1]  &&  any( xi < TF$domain[1, ] ) ) next

        if( domaincheck[2]  &&  any( TF$domain[2, ] < xi ) ) next

        out[i, ] =  TF$fun( xi )
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



validate.elemTransferFunction  <-  function( TF, points=1300, tol=5.e-7, domain=NULL )
    {
    out     = TRUE
    mess    = character(0)

    ok  = is.numeric(points)  &&  length(points)==1  &&  8<=points
    if( ! ok )
        {
        log_string( ERROR, "points='%s' is invalid.", as.character(points)[1] )
        return(NULL)
        }

    m   = dimension(TF)

    if( is.null(domain) )
        {
        domain  = TF$domain
        }
    else
        {
        if( ! is.matrix(domain) )
            {
            dim(domain) = NULL
            domain      = prepareNxM( domain, 2 )
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



    if( m == 1 )
        {
        interval    = as.numeric( domain )
        finite      = is.finite(interval)

        if( all( finite == c(TRUE,FALSE) ) )
            interval[2] = interval[1] + 10
        else if( all( finite == c(FALSE,TRUE) ) )
            interval[1] = interval[2] - 10
        else if( all( finite == c(FALSE,FALSE) ) )
            interval    = c(-10,10)

        #   print( interval )

        x   = seq( interval[1], interval[2], len=points )
        y   = TF$fun( x )

        #   check monotonicity
        delta   = diff(y)
        ok  = all(delta<0)  ||  all(0<delta)
        if( ! ok )
            {
            mess = c( mess, "The function is not monotone." )
            out = FALSE
            }

        # check that all y are inside the given range
        distance    = pmax( TF$range[1]-y, y-TF$range[2], 0 )
        if( any(0 < distance) )
            {
            idx = which.max( distance )
            mess = c( mess, sprintf( "The function maps %d of %d points out of the range, e.g. f(%g)=%g which is not in [%g,%g].  dist=%g.",
                                    sum(0<distance), length(distance), x[idx], y[idx], TF$range[1], TF$range[2], distance[idx] ) )
            out = FALSE
            }

        line = sprintf( "range-test points = %d, max(distance)=%g.\n", length(distance), max(distance)  )
        cat( line )

        #   check that NA returns NA
        z = TF$fun( NA_real_ )
        if( ! is.na(z) )
            {
            mess = c( mess, sprintf( "The function does not handle NA correctly; it returns %g instead of NA.", z ) )
            out = FALSE
            }

        if( is.invertible(TF) )
            {
            x.back  = TF$funinv(y)
            delta   = abs( x.back - x )   #; print( max(delta) )

            side    = interval[2] - interval[1]
            invalid = (side * tol  <  delta)

            if( any(invalid,na.rm=TRUE) )
                {
                mess = c( mess, sprintf( "The function failed the round-trip invertibility tolerance test, for %d of %d points.", sum(invalid,na.rm=TRUE), length(invalid) ) )
                mess = c( mess, sprintf( "    max(deltarow)=%g > %g = %g * %g = side * tol.", max(delta,na.rm=TRUE), side*tol, side, tol ) )
                idx = which.max(delta)
                x0  = x[idx]
                y0  = y[idx]
                x0.back = x.back[idx]
                mess = c( mess, sprintf( "    %g  ->  %g  ->  %g", x0, y0, x0.back ) )
                out = FALSE
                }

            mess = c( mess, sprintf( "round-trip points = %d, max(delta)=%g.\n", length(x), max(delta)  ) )
            }
        }
    else
        {
        #mess    = sprintf( "dimension %d validation not working.", dimension(TF) )
        #   how many points on side of the cube
        n   = pmax( ceiling( points^(1/m) ), 2 )

        grid    = vector( m, mode='list' )
        for( k in 1:m )
            grid[[k]]   = seq( domain[1,k], domain[2,k], length.out=n )

        x           = as.matrix( expand.grid(grid) )  #; print( str(x) )
        colnames(x) = colnames(domain)
        y           = x
        valid       = logical(nrow(x))
        #   disti       = numeric(m)
        distance    = numeric(nrow(x))

        for( i in 1:nrow(x) )
            {
            y[i, ]      = TF$fun( x[i, ] )

            disti   = pmax( TF$range[1, ]-y[i, ], y[i, ]-TF$range[2, ], 0 )

            distance[i] = max( disti )
            }

        if( ! all( is.finite(y) ) )
            {
            invalid = ! is.finite( .rowSums( y, nrow(y), ncol(y) ) )

            mess = c( mess, sprintf( "The transfer function maps %d of %d test inputs to a non-finite value.", sum(invalid), length(invalid) ) )
            out = FALSE
            }

        if( any(0 < distance) )
            {
            idx = which.max( distance )
            mess = c( mess, sprintf( "The function maps %d of %d points outside the range.  maxdist=%g.",
                                    sum(0<distance), length(distance), distance[idx] ) )
            x0  = x[idx, ]
            y0  = y[idx, ]
            mess = c( mess, sprintf( "    %s  ->  %s", nicevector(x0), nicevector(y0) ) )
            out = FALSE
            }

        mess = c( mess, sprintf( "range-test points = %d, max(distance to range)=%g.\n", length(distance), max(distance)  )  )


        if( is.invertible(TF) )
            {
            x.back  = y
            for( i in 1:nrow(x) )
                x.back[i, ]  = TF$funinv( y[i, ] )

            if( ! all( is.finite(x.back) ) )
                {
                invalid = ! is.finite( .rowSums( x.back, nrow(x.back), ncol(x.back) ) )
                mess = c( mess, sprintf( "The transfer function inverse maps %d of %d test inputs to a non-finite value.", sum(invalid), length(invalid) ) )
                idx = which( invalid )[1]
                x0  = x[idx, ]
                y0  = y[idx, ]
                x0.back = x.back[idx, ]
                mess = c( mess, sprintf( "    %s  ->  %s  ->  %s", nicevector(x0), nicevector(y0), nicevector(x0.back) ) )
                out = FALSE
                }

            delta   = abs( x.back - x )     #; print( max(delta) )

            #   relativize the columns of delta
            #for( k in 1:m )
            #    delta[ ,k]  = delta[ ,k] / (domain[2,k] - domain[1,k])

            #print( mess )
            #print( max(delta) )

            sidemax     = max( domain[2, ] - domain[1, ] )

            deltarow    = apply( delta, 1, max )        #.rowSums( delta, nrow(delta), ncol(delta), na.rm=TRUE ) / m

            invalid = (sidemax * tol  <  deltarow)

            if( any(invalid,na.rm=TRUE) )
                {
                mess = c( mess, sprintf( "The function failed the round-trip invertibility tolerance test, for %d of %d points.", sum(invalid,na.rm=TRUE), length(invalid) ) )
                mess = c( mess, sprintf( "    max(deltarow)=%g > %g = %g * %g = sidemax * tol.", max(deltarow,na.rm=TRUE), sidemax*tol, sidemax, tol ) )
                idx = which.max(deltarow)
                x0  = x[idx, ]
                y0  = y[idx, ]
                x0.back = x.back[idx, ]
                mess = c( mess, sprintf( "    %s  ->  %s  ->  %s", nicevector(x0), nicevector(y0), nicevector(x0.back) ) )
                out = FALSE
                }

            mess = c( mess, sprintf( "round-trip points = %d, max(deltarow)=%g.\n", nrow(x), max(deltarow,na.rm=TRUE) )  )
            }

        # out     = FALSE
        }

    if( 0 < length(mess) )  attr( out, 'message' ) = mess

    return(out)
    }




print.elemTransferFunction  <-  function( x, id=NULL, ... )
    {
    #   print( sys.status() )

    #for( where in -length(sys.frames()):length(sys.frames()) )
    #    {
    #    test = deparse( substitute(x,sys.frame(where)) )
    #    print(test)
    #    }


    # theName = deparse( substitute(x) )   # theName is 'x' when print() is called through UseMethod(),  hmm....    #,parent.frame(1)

    if( is.null(id) )
        id = "This"
    else
        id = id[1]

    n   = dimension(x)
    if( n == 1 )
        {
        cat( sprintf( "%s is a univariate TransferFunction.\n", id  ) )
        cat( sprintf( "domain:      [%g,%g]  (%s)\n", x$domain[1], x$domain[2], colnames(x$domain) ) )
        cat( sprintf( "range:       [%g,%g]  (%s)\n", x$range[1], x$range[2], colnames(x$range) ) )
        }
    else
        {
        cat( sprintf( "%s is a multivariate TransferFunction, of dimension %d.\n", id, n ) )
        cat( "domain:\n" )
        print( x$domain )
        cat( "range:\n" )
        print( x$range )
        }

    cat( sprintf( "invertible:  %s\n", ifelse( is.invertible(x), 'Yes', 'No' ) ) )

    oclass  = orientation(x)
    if( is.na(oclass) || oclass==0 )
        oclass  = 'NA'
    else
        oclass  = ifelse( 0<oclass, 'preserving', 'reversing' )
    cat( sprintf( "orientation: %s\n", oclass ) )


    pass    = validate( x )
    cat( sprintf( "validation:  %s\n", ifelse( pass, 'Passed', 'Failed' ) ) )
    if( ! pass )
        cat( attr(pass,'message'), sep='\n' )

    return( invisible(TRUE) )
    }



orientation.elemTransferFunction <- function(TF)
    {
    n   = dimension(TF)

    theDomain   = TF$domain

    x   = matrix( theDomain[1, ], n+1, n, byrow=TRUE )

    for( k in 1:n )
        x[k,k]    = theDomain[2,k]

    #   print(x)

    y   = transfer( TF, x )                                 #; print(y)
    if( is.null(y) )    return(NA_real_)

    mat = y[ 1:n, ] - matrix( y[n+1, ], n, n, byrow=TRUE )  #; print(mat)

    return( det(mat) )
    }


nicevector <- function( v )
    {
    paste0( '(', paste0( sprintf( "%g", v ), collapse=',' ), ')', collapse='' )
    }



#   TF      a transfer function that has already been validated
makeInverseTF  <-  function( TF, domain=c(0,1) )
    {
    #   estimate the rough gamma of TF
    range   = TF( domain )

    xmid    = 0.5*sum(domain)
    ymid    = TF( xmid )
    gamma   = log( (ymid - range[1])/(range[2] - range[1]) ) / log(0.5)  #; print(gamma)

    x   = seq( 0, 1, by=1/512 ) ^ (1/gamma)
    x   = (1-x)*domain[1]  +  x*domain[2]
    y   = TF(x)

    out = splinefun( y, x, method='monoH.FC' )

    return( out )
    }
