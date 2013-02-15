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

float integ_sumsquares_r(int,float**,int*,float,float);
float integ_sumsquares_z(int,int,float**,int*,float);
void get_diff_elems_and_volume(int*,float&,float&,float&,float&);

//============================================================================//
//============================================================================//
float integ_sumsquares_3D(float **vec, int *dims, bool volume_avg) {
    int k;
    float integ, vec2_k1, vec2_k2, dz, dr, dp, V;
    get_diff_elems_and_volume(dims,dz,dr,dp,V);
    integ = 0.0;
    for(k=1; k<dims[2]; k++) {
        vec2_k1 = integ_sumsquares_r(k-1,vec,dims,dr,dz);
        vec2_k2 = integ_sumsquares_r(k  ,vec,dims,dr,dz);
        integ += (vec2_k1+vec2_k2)/2. * dp;
    }
    vec2_k1 = integ_sumsquares_r(dims[2]-1,vec,dims,dr,dz);
    vec2_k2 = integ_sumsquares_r(0        ,vec,dims,dr,dz);
    integ += (vec2_k1+vec2_k2)/2. * dp;
    if(volume_avg)
        return integ/V;
    else
        return integ;
}

//============================================================================//
float integ_sumsquares_2D(float **vec, int *dims, bool volume_avg) {
    float integ, dz, dr, dp, V;
    get_diff_elems_and_volume(dims,dz,dr,dp,V);
    integ = 2.*pi*integ_sumsquares_r(0,vec,dims,dr,dz);
    if(volume_avg)
        return integ/V;
    else
        return integ;
}

//============================================================================//
float integ_sumsquares_r(int k, float **vec, int *dims, float dr, float dz) {
    float integ, vec2_j1, vec2_j2;
    integ = 0.0;    
    for(int j=1; j<dims[1]; j++) {
        vec2_j1 = integ_sumsquares_z(k,j-1,vec,dims,dz);
        vec2_j2 = integ_sumsquares_z(k,j  ,vec,dims,dz);
        integ += ((vec2_j1+vec2_j2)/2. * dr) * ((2.*j-1.)/2. * dr); // r*dr
    }
    return integ;
}

//============================================================================//
float integ_sumsquares_z(int k, int j, float **vec, int *dims, float dz) {
    int n1, n2;
    float integ, vec2_i1, vec2_i2;
    integ = 0.0;
    for(int i=1; i<dims[0]; i++) {
        n1 = fn(i-1,j,k,dims[0],dims[1]);
        n2 = fn(i  ,j,k,dims[0],dims[1]);
        vec2_i1 = 0.0;
        vec2_i2 = 0.0;
        for(int m=0; m<ndims; m++) {
            vec2_i1 += vec[m][n1]*vec[m][n1];
            vec2_i2 += vec[m][n2]*vec[m][n2];
        }
        integ += (vec2_i1+vec2_i2)/2. * dz;
    }    
    return integ;
}

//============================================================================//
void get_diff_elems_and_volume(int *dims, float &dz, float &dr, float &dp, 
                                float &V) {
    dz = 2.*42.2222/(dims[0]-1);
    dr = 28.1944/(dims[1]-1);
    dp = 2*pi/dims[2];
    V  = pi*(28.1944*28.1944)*(2.*42.2222);
}

//============================================================================//
//============================================================================//
void Fourier_Decomp(float **vec, int *dims, int i, int j, int vec_comp,
                    float &c0, float &c1) {
    // Compute n=0 and n=1 Fourier decomposition of a cylindrical vector
    // For c1, k0 is a scan over the various starting positions for theta0, 
    // which is an arbitrary choice.  Thus the various theta0 results are 
    // averaged to give an unbiased decomposition.
    int k, k0, n;
    float theta_k, c1_r, c1_i, c1_sum;
    
    c0 = 0.0;
    for(k=0; k<dims[2]; k++) {
        n = fn(i,j,k,dims[0],dims[1]);
        c0 += vec[vec_comp][n];
    }
    c0 = sqrt(c0*c0)/(1.*dims[2]);    
   
    c1_sum = 0.0;
    for(k0=0; k0<dims[2]; k0++) {
        c1_r = 0.0;
        c1_i = 0.0;
        for(k=0; k<dims[2]; k++) {
            theta_k = 2*pi*(k-k0)/(1.*dims[2]);
            n = fn(i,j,k,dims[0],dims[1]);
            c1_r += vec[vec_comp][n]*cos(theta_k);
            c1_i += vec[vec_comp][n]*sin(theta_k);
        }
        c1_sum += sqrt(2.*c1_r*c1_r + 2.*c1_i*c1_i)/(1.*dims[2]);
    }
    c1 = c1_sum/(1.*dims[2]);
}

//============================================================================//
//============================================================================//


