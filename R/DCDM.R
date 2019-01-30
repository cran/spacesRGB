


#   SMPTE 431-2 defines the DCDM color encoding equations. 
#   The equations for the decoding of the encoded color information are the 
#   inverse of the encoding equations
#   Note: Here the 4095 12-bit scalar is not used since the output of CTL is 0-1.

#   XYZ     -> X'Y'Z'    
dcdm_encode <- function( XYZ )
    {
    XYZp    = ( (48./52.37) * XYZ ) ^ (1/2.6)

    return( XYZp )
    }



#   X'Y'Z'    -> XYZ
dcdm_decode <- function( XYZp )
    {
    XYZ =   (52.37/48.0) * XYZp ^ 2.6

    return( XYZ )
    }

