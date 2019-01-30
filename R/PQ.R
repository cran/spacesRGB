
#   For the equations here, see:
#   https://en.wikipedia.org/wiki/High-dynamic-range_video#Perceptual_Quantizer

#   Lmax    maximum luminance, nit
#
#   There is one modification, the input display signal V is increased (by about 1.e-6) if necessary so that fun(V) is defined (and 0).
#   See the call to pmax() near the bottom.

PQ.EOTF  <-  function( Lmax=10000 )
    {
    ok  = is.numeric(Lmax)  &&  length(Lmax)==1  &&  0<Lmax
    if( ! ok )
        {
        log.string( ERROR, "Lmax='%s' is invalid.", as.character(Lmax) )
        return(NULL)
        }   
        
    #   These d's are the c's multiplied by 128
    #   This rescaling is probably unnecessary because 128 is a power of 2, so the c's are dyadic rationals.
    #   But do this anyway because it looks so much cleaner
    d1  =  107          # const float pq_c1 = 0.8359375; // 3424.0 / 4096.0 or pq_c3 - pq_c2 + 1.0;
    d2  = 2413          # const float pq_c2 = 18.8515625; // ( 2413.0 / 4096.0 ) * 32.0;
    d3  = 2392          # const float pq_c3 = 18.6875; // ( 2392.0 / 4096.0 ) * 32.0;    
    d4  =  128          # pq_c4 is 1
        
    m1  = 652.5 / 4096  # const float pq_m1 = 0.1593017578125; // ( 2610.0 / 4096.0 ) / 4.0;
    m2  = 2523 / 32     # const float pq_m2 = 78.84375; // ( 2523.0 / 4096.0 ) * 128.0;
    
    #   Vmin    = (d1/d4)^m2
    Vmin    = 0
    
    domain  = matrix( c(Vmin,1), 2, 1, dimnames=list(NULL,"non-linear signal") )
    
    #   make range extra small so composition works
    range   = matrix( c(0,Lmax), 2, 1, dimnames=list(NULL,"linear display") )
    
    fun     <-  function(V) { Vm = V^(1/m2) ; Lmax * ( pmax(d4*Vm - d1,0)/(-d3*Vm + d2) )^(1/m1) }      # note the pmax() here !
    funinv  <-  function(L) { L0 = (L/Lmax)^m1 ; ((d2*L0 + d1)/(d3*L0 + d4))^m2 }
    
    TransferFunction( fun, funinv, domain, range, id=sigfunction() )
    }
        