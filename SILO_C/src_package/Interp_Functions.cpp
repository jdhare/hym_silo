//============================================================================//
/*

Clayton Myers
Interp_Functions.cpp
Created:  28 July 2010

Interpolation support for 1D and 2D arrays

*/
//============================================================================//
//============================================================================//

#include <HYM_SILO.hpp>
void find_nearest(float*,int,float,int&,int&);

//============================================================================//
float interp_1d(float *xg, int Nx, float *yg, float x) {
    int iL, iR;

    find_nearest(xg,Nx,x,iL,iR);
    
    float x0 = xg[iL];
    float x1 = xg[iR];
    float y0 = yg[iL];
    float y1 = yg[iR];
    
    return y0 + (x-x0)*(y1-y0)/(x1-x0);
}

//============================================================================//
float interp_1d_pts(float x, float x0, float x1, float y0, float y1) {
    return y0 + (x-x0)*(y1-y0)/(x1-x0);
}

//============================================================================//

float interp_2d(float *zg, int Nz, float *rg, int Nr, float *f, 
                float z, float r) {
    int iL, iR, jL, jR;
    float x, y, x1, x2, y1, y2, f_Q11, f_Q12, f_Q21, f_Q22, f_R1, f_R2; 

    find_nearest(zg,Nz,z,iL,iR);
    find_nearest(rg,Nr,r,jL,jR);

    x = z;
    y = r; 
    x1 = zg[iL];
    x2 = zg[iR];
    y1 = rg[jL];
    y2 = rg[jR];
    
    f_Q11 = f[iL+jL*Nz];
    f_Q12 = f[iL+jR*Nz];
    f_Q21 = f[iR+jL*Nz];
    f_Q22 = f[iR+jR*Nz];
    
    f_R1 = ((x2-x)*f_Q11 + (x-x1)*f_Q21)/(x2-x1);
    f_R2 = ((x2-x)*f_Q12 + (x-x1)*f_Q22)/(x2-x1);
    return ((y2-y)*f_R1  + (y-y1)*f_R2 )/(y2-y1);
}
    
//============================================================================//

void find_nearest(float *xg, int Nx, float x, int &iL, int &iR) {
    // Find the nearest indices (left and right) using the bisection method.
    int iMid;
    iL = 0;
    iR = Nx-1;
    while((iR-iL) > 1) {
        iMid = (int)(round((iR+iL)/2.));
        if(xg[iMid] < x)
            iL = iMid;
        else   
            iR = iMid;
    }
}

//============================================================================//
//============================================================================//

