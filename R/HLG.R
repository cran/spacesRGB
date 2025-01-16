
#   OOTF            from scene linear to display linear
#   Lb and Lw       interval for the display linear

HLG.OOTF  <-  function( gamma=1.2, Lb=0, Lw=1000 )
    {
    ok  = is.numeric(gamma)  &&  length(gamma)==1  &&  1<=gamma    
    if( ! ok )
        {
        log_level( ERROR, "gamma = '%s' is invalid.", as.character(gamma) )
        return(NULL)
        }       
    
    ok  = length(Lb)==1  &&  length(Lw)==1  &&  (0 <= Lb)  &&  (Lb < Lw)
    if( ! ok )
        {
        log_level( ERROR, "Lb=%g  and/or  Lw=%g are invalid.", Lb, Lw )
        return(NULL)
        }    
    
    alpha   = Lw - Lb
    
    weights.rgb = c(0.2627,0.6780,0.0593)
    
    # scene linear to display linear
    fun <- function( rgb )  
        {
        Ys  = sum( weights.rgb * rgb )  # scene luminance, in [0,1]
        
        rgb = Ys^(gamma-1) * rgb   #; print(rgb)   # rgb is in [0,1]^3
        
        return( rgb*Lw  +  (1-rgb)*Lb )
        }
    
    # display linear to scene linear
    funinv <- function( rgb )
        {
        Yd  = sum( weights.rgb * rgb )
        
        s   = ((Yd - Lb)/alpha)
        
        if( s <= 0 )    return( numeric(3) )    # all 0s
        
        s   = s ^ ((1-gamma)/gamma) 
        
        return( s * (rgb - Lb) / alpha )
        }
    
    cnames  = paste( "linearscene", c('.R','.G','.B'), sep='' )
    domain  = matrix( c(0,1),   2, 3, dimnames=list(NULL,cnames) )
    
    cnames  = paste( "lineardisplay", c('.R','.G','.B'), sep='' )
    range   = matrix( c(Lb,Lw), 2, 3, dimnames=list(NULL,cnames) )

    TransferFunction( fun, funinv, domain, range, id=sigfunction() )
    }
    
    
    
#   OETF    from scene linear to display signal
    
    
roota <- function()
    {
    fun <- function( x )    { x * log( (x + 11/4)/x ) - 1/2 }
    
    uniroot( fun, c(0.17,0.18), tol = .Machine$double.eps^0.75 )$root
    }
    
    
HLG.OETF <- function()
    {
    a   = roota()      #; print(a)  # 0.178832772656984     # from roota() above.  This takes only ~80 usec.     # 0.17883277
    b   = 1 - 4*a
    c   = 0.5 - a*log(4*a)
    
    fun <- function( lin )
        {
        lin = as.numeric(lin)
        
        sig = lin   # to get the size right, and any NAs in lin[] get copied over too :-)
        
        lo  = lin < 1/12    
        lo[ is.na(lo) ] = FALSE     #; print(lo)
        sig[lo] = sqrt( 3*lin[lo] )
        
        hi  = ! lo  #   &  lin < 1 
        hi[ is.na(hi) ] = FALSE 
        sig[hi] = a*log(12*lin[hi] - b) + c
        
        end = lin == 1
        end[ is.na(end) ] = FALSE
        sig[end] =  1
        
        return(sig)
        }
    
    funinv  <- function( sig )
        {
        sig = as.numeric(sig)
        
        lin = sig   # to get the size right
        
        lo  = sig <= 1/2
        lo[ is.na(lo) ] = FALSE        
        lin[lo] = sig[lo]^2 / 3
        
        hi  = ! lo  #  &  sig < 1 
        hi[ is.na(hi) ] = FALSE        
        lin[hi] = (exp((sig[hi]-c)/a) + b) / 12
        
        end = sig == 1
        end[ is.na(end) ] = FALSE        
        lin[end]=  1
        
        return(lin)
        }
    
    linmax  = 1
    sigmax  = fun(linmax)
    
    domain  = matrix( c(0,linmax), 2, 1, dimnames=list(NULL,"scene linear") )
    range   = matrix( c(0,sigmax), 2, 1, dimnames=list(NULL,"display signal") )
    
    TransferFunction( fun, funinv, domain, range, id='HLG.OETF' )
    }
        

    

    
    
    
#   HLG2D.OOTF  is just for testing and plotting a 2D example
#
#   weights     should sum to 1
#
HLG2D.OOTF  <-  function( gamma=1.2, weights=c(0.5,0.5) )
    {
    ok  = is.numeric(gamma)  &&  length(gamma)==1  &&  1<=gamma    
    if( ! ok )
        {
        log_level( ERROR, "gamma = '%s' is invalid.", as.character(gamma) )
        return(NULL)
        }       
        
    ok  = is.numeric(weights)  &&  length(weights)==2  &&  all(0 < weights)  &&  abs(sum(weights) - 1) < 1.e-12
    if( ! ok )
        {
        log_level( ERROR, "weights = '%s' is invalid.", paste0( weights,collapse=',') )
        return(NULL)
        } 
        
        
    #   alpha   = 1

    # scene linear to display linear
    fun <- function( rg )  
        {
        Ys  = sum( weights * rg )
        
        s   = Ys^(gamma-1)      # s is in [0,1]
        
        return( s*rg )
        }
    
    # display linear to scene linear
    funinv <- function( rg )
        {
        Yd  = sum( weights * rg )
        
        s   = Yd
        
        if( s <= 0 )    return( numeric(2) )    # all 0s
        
        s   = s ^ ((1-gamma)/gamma) 
        
        return( s * rg )
        }
    
    cnames  = paste( "linearscene", c('.R','.G'), sep='' )
    domain  = matrix( c(0,1),   2, 2, dimnames=list(NULL,cnames) )
    
    cnames  = paste( "lineardisplay", c('.R','.G'), sep='' )
    range   = matrix( c(0,1), 2, 2, dimnames=list(NULL,cnames) )

    TransferFunction( fun, funinv, domain, range, id=sigfunction() )
    }
        