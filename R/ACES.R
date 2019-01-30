
####----    Chromaticities of some common primary sets    ----####


AP0_PRI     = NULL      #   ACES Primaries from SMPTE ST2065-1
AP1_PRI     = NULL      #   Working space and rendering primaries for ACES 1.0
REC709_PRI  = NULL
REC2020_PRI = NULL
P3D65_PRI   = NULL
P3D60_PRI   = NULL
P3DCI_PRI   = NULL

p.AP0_2_XYZ_MAT = NULL
p.XYZ_2_AP0_MAT = NULL
p.AP1_2_XYZ_MAT = NULL
p.XYZ_2_AP1_MAT = NULL

p.AP0_2_AP1_MAT = NULL
p.AP1_2_AP0_MAT = NULL
p.AP1_RGB2Y     = NULL

TINY    = 1e-10


#   this should be called from .onLoad()
makeAllPrimaries <- function()
    {
    #   ACES Primaries from SMPTE ST2065-1
    AP0_PRI     <<- makePrimaries4x2( c(0.73470,0.26530,       0,      1, 0.0001,-0.07700, 0.32168,0.33767 ) )

    #   Working space and rendering primaries for ACES 1.0
    AP1_PRI     <<- makePrimaries4x2( c( 0.71300,0.29300, 0.16500,0.83000, 0.12800,0.04400, 0.32168,0.33767 ) )

    REC709_PRI  <<- makePrimaries4x2( c( 0.64000,0.33000, 0.30000,0.60000, 0.15000,0.06000, 0.31270,0.32900 ) )

    REC2020_PRI <<- makePrimaries4x2( c( 0.70800,0.29200, 0.17000,0.79700, 0.13100,0.04600, 0.31270,0.32900 ) )

    P3D65_PRI   <<- makePrimaries4x2( c( 0.68000,0.32000, 0.26500,0.69000, 0.15000,0.06000, 0.31270,0.32900 ) )

    P3D60_PRI   <<- makePrimaries4x2( c( 0.68000,0.32000, 0.26500,0.69000, 0.15000,0.06000, 0.32168,0.33767 ) )

    P3DCI_PRI   <<- makePrimaries4x2( c( 0.68000,0.32000, 0.26500,0.69000, 0.15000,0.06000, 0.31400,0.35100 ) )
        
    return(TRUE)
    }


makePrimaries4x2 <- function( xy )
    {
    out = base::matrix( xy, 4, 2, byrow=TRUE )
    
    rownames(out)   = c('R','G','B','W')
    colnames(out)   = c('x','y')
    
    return(out)
    }




rgb_2_saturation    <- function( rgb )
    {
    #   note that multiplying rgb by a constant does not change saturation, except very near 0
    return( ( max(rgb,TINY) - max(min(rgb),TINY)) / max(rgb,1e-2) )
    }

rgb_2_hue   <- function( rgb )
    {
    #   Returns a geometric hue angle in degrees (0-360) based on RGB values.
    #   For neutral colors, hue is undefined and the function will return a quiet NaN value.

    if (rgb[1]==rgb[2]  &&  rgb[2]==rgb[3] )   
        #   RGB triplets where RGB are equal have an undefined hue
        return( 0 )   #NaN )

    hue = (180/pi) * atan2( sqrt(3)*(rgb[2]-rgb[3]), 2*rgb[1]-rgb[2]-rgb[3] )

    if( hue < 0 ) hue = hue + 360

    return( hue )
    }
    
    
rgb_2_yc    <- function( rgb, ycRadiusWeight=1.75 )
    {
    #   Converts RGB to a luminance proxy, here called YC
    #   YC is ~ Y + K * Chroma
    #   Constant YC is a cone-shaped surface in RGB space, with the tip on the 
    #   neutral axis, towards white.
    #   YC is normalized: RGB 1 1 1 maps to YC = 1
    #   
    #   ycRadiusWeight defaults to 1.75, although can be overridden in function 
    #   call to rgb_2_yc
    #   ycRadiusWeight = 1 -> YC for pure cyan, magenta, yellow == YC for neutral 
    #   of same value
    #   ycRadiusWeight = 2 -> YC for pure red, green, blue  == YC for  neutral of 
    #   same value.

    r = rgb[1]
    g = rgb[2]
    b = rgb[3]
  
    chroma = sqrt(b*(b-g)+g*(g-r)+r*(r-b))

    return( ( b + g + r + ycRadiusWeight * chroma) / 3 )
    }
    
    
calc_sat_adjust_matrix  <- function( sat, rgb2Y )
    {
    #   This function determines the terms for a 3x3 saturation matrix that is
    #   based on the luminance of the input.

    M   = (1-sat) * matrix( rgb2Y, 3, 3, byrow=TRUE )  +  sat * diag(3)

    # M = transpose_f33(M);   Do *not* transpose, because this R package uses the standard matrix convention
  
    return( M )
    }


calc_sat_adjust_matrix_long_way  <- function( sat, rgb2Y )
    {
    #   This function determines the terms for a 3x3 saturation matrix that is
    #   based on the luminance of the input.

    M   = matrix( 0, 3, 3 )
    
    M[1,1] = (1 - sat) * rgb2Y[1] + sat
    M[2,1] = (1 - sat) * rgb2Y[1]
    M[3,1] = (1 - sat) * rgb2Y[1]

    M[1,2] = (1 - sat) * rgb2Y[2]
    M[2,2] = (1 - sat) * rgb2Y[2] + sat
    M[3,2] = (1 - sat) * rgb2Y[2]

    M[1,3] = (1 - sat) * rgb2Y[3]
    M[2,3] = (1 - sat) * rgb2Y[3]
    M[3,3] = (1 - sat) * rgb2Y[3] + sat

    # M = transpose_f33(M);   Do *not* transpose, because this R package uses the standard matrix convention
  
    return( M )
    }

    
    
#---- Functions to compress highlights ----#
# allow for simulated white points without clipping

#   Y           vector of color values to adjust (white scaled to around 1.0)
#   new_wht     white adjustment (e.g. 0.9 for 10% darkening)
#   width       adjusted width (e.g. 0.25 for top quarter of the tone scale)

roll_white_fwd <- function( Y, new_wht, width )
    {
    x0  = -1.0
    x1  = x0 + width
    y0  = -new_wht
    y1  = x1
    m1  = (x1 - x0)
    a   = y0 - y1 + m1
    b   = 2 * ( y1 - y0) - m1
    #   c   = y0
    t   = (-Y - x0) / (x1 - x0)     # t is now a vector
    
    out = Y
    
    mask0       = t < 0
    out[mask0]  = -(t[mask0] * b + y0)

    mask1       = !mask0  &  t <= 1
    out[mask1]  = -(( t[mask1] * a + b) * t[mask1] + y0)

    #   and if t > 1 then out =  Y
    
    return( out )
    }    


    
#   Y           vector of color values to adjust (white scaled to around 1.0)
#   new_wht     white adjustment (e.g. 0.9 for 10% darkening)
#   width       adjusted width (e.g. 0.25 for top quarter of the tone scale)

roll_white_rev <- function( Y, new_wht, width )
    {
    x0  = -1
    x1  = x0 + width
    y0  = -new_wht
    y1  = x1
    m1  = (x1 - x0)
    a   = y0 - y1 + m1
    b   = 2 * ( y1 - y0) - m1
    # float c = y0
    
    out = Y
    
    mask0   = -Y < y0
    
    out[mask0] = -x0[mask0]
    
    mask1   = !mask0  &  -Y <= y1
    
    if( any(mask1) )
        {
        c           = y0 + Y[mask1]
        discrim     = sqrt( b * b - 4. * a * c)
        t           = ( 2. * c) / ( -discrim - b)
        out[mask1]  = -(( t * ( x1 - x0)) + x0)
        }

    return( out )
    }
    