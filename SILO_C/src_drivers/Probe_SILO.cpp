//============================================================================//
/*

Clayton Myers
Probe_SILO.cpp
Created:  18 August 2009
Modified: 30 November 2009

Main function file for reading magnetics data for simulated SSX probe
diagnostics from HYM's .silo storage format.  This program calculates the index
coordinates in HYM that correspond to the physical locations of the magnetic
probes in SSX.
    
The probe locations are given fractionally:
    In Z: (FC Length =  12.0 in = 61.0 cm)
        East:     -21.6 cm  ->  8.9/61.0 = 0.146
        Midplane:   0.0 cm  -> 30.5/61.0 = 0.500
        West:     +21.6 cm  -> 52.1/61.0 = 0.854
    
    In R: (FC Radius = 8.0 in = 20.3 cm)
        Probe 1:    1.0 in  ->   1.0/8.0 = 0.125
        Probe 2:    2.0 in  ->   2.0/8.0 = 0.250
        Probe 3:    3.0 in  ->   3.0/8.0 = 0.375
                ......................
        Probe 8:    8.0 in  ->   8.0/8.0 = 1.000

    In Phi: (Could be arbitrarily rotated)
        Angle 1:   45.0 deg ->    45/360 = 0.125
        Angle 2:  135.0 deg ->   135/360 = 0.375
        Angle 3:  225.0 deg ->   225/360 = 0.625
        Angle 4:  315.0 deg ->   315/360 = 0.875
        
The location arrays (e.g. zLoc) contain lists of the fractional locations of the
probes along the relevant dimension (see above).
        
The index arrays (e.g. zInd) contain lists of the calculated nearest indices
to the probe locations in the mesh along the relevant dimension.

The interpolation arrays (e.g. zInterp) contain lists of the fractional residual
distance between the actual location of the probe and the integer node to which 
it was mapped.  Thus, these values have a range of [-0.5,+0.5].  These values 
can be used to interpolate among several nodes to get more accurate data values
at the probe locations if desired.  A residual of zero (0.000) means that there 
is an exact mapping between the probe location and the mesh.

The radial dimension (r) requires divisibility by 8 to be residual-free and the 
phi dimension (phi) requires divisibility by 4 to be residual-free.  These two 
criteria are often satisfied, so only the axial (z) values need to be 
interpolated in most cases.

Dimensions of HYM mesh in (z,r,phi) format:
  int dims[ndims] = {129,65,32};  // CoH Mesh
  int dims[ndims] = {257,129,16}; // CounterH Mesh
  
----------------------------------------

Command line arguments:
    (1) silo_path (string) -- path to the location of the silo database from 
            which the probe data is being extracted.
    (2) probe_path (string) -- path to the location of the output file for the
            extracted probe data.
    (3) cycle (int) -- cycle number for extracted data.  If cycle == 0, data 
            will be extracted from all available cycles.
    (4) phi_rot (float) -- Arbitrary rotation angle (in degrees) for the 
            extracted probe array.
            
*/
//============================================================================//
//============================================================================//

#include <HYM_SILO.hpp>
#include <SILO_Read.hpp>

//============================================================================//
//============================================================================//
void ReadArgs(int,char**,char*&,char*&,int&,float&);
void Write_ProbeData(int,char*,char*);
void Map_Indices(float*,int,int,int*,float*);
void Print_Indices(int,char*,int*,float*);
void Reset_pLoc_Interval();

// Number of probes in each dimension
const int zLen = 3;
const int rLen = 8;
const int pLen = 4;

// Fractional probe locations
float zLoc[zLen] = {0.146,0.500,0.854};
// float zLoc[zLen] = {0.235,0.500,0.765}; // For shifted probe locations
float rLoc[rLen] = {0.125,0.250,0.375,0.500,0.625,0.750,0.875,1.000};
float pLoc[pLen] = {0.00,0.25,0.50,0.75};

char *silo_name = "HYM";
char *stopmsg = "Stopping probe data extraction.";

//============================================================================//
int main(int argc, char *argv[]) {
    int cycle, Ncyc;
    float phi_rot;
    char *silo_path=NULL, *probe_path=NULL;
    
    // Read the command line arguments:
    ReadArgs(argc,argv,silo_path,probe_path,cycle,phi_rot);
    
    // Rotate the phi locations by phi_rot:
    for(int k=0; k<pLen; k++)
        pLoc[k] += phi_rot/360.0;
    Reset_pLoc_Interval();
    
    // Determine the available number of cycles:
    // Ncyc = Get_Ncyc(silo_path,stopmsg);
    Ncyc = 250;

    // Write (to a file) the probe data for the requested cycle argument:
    if(cycle == 0) {
        for(int m=1; m <= Ncyc; m++)
            Write_ProbeData(m,silo_path,probe_path);
    }
    else if(cycle >= 1 && cycle <= Ncyc)
        Write_ProbeData(cycle,silo_path,probe_path);
    else {            
        char message[1001];
        sprintf(message,"  %s\n      %s%d%s\n      %s%d%s\n  %s",
                "Error in command line arguments:",
                "The requested cycle number is not valid (cycle = ",cycle,").",
                "The cycle number must be in the range [1,",Ncyc,"].",
                stopmsg);
        StopExecution(message);
    }
}

//============================================================================//
//============================================================================//
void ReadArgs(int argc, char **argv, char *&silo_path, char *&probe_path,
              int &cycle, float &phi_rot) {
    // Count the initial command line arguments:
    if(!(argc == (1+4))) {
        char message[1001];
        sprintf(message,"  %s\n  %s",
                "An improper number of command line arguments was found.",
                stopmsg);
        StopExecution(message);
    }
        
    // Distribute the command line arguments
    silo_path = argv[1];
    probe_path = argv[2];
    VerifyPath(silo_path,stopmsg);
    VerifyPath(probe_path,stopmsg);
    ConvertToInt(argv[3],cycle,stopmsg);
    ConvertToFloat(argv[4],phi_rot,stopmsg);
}

//============================================================================//
void Write_ProbeData(int cycle, char *silo_path, char *probe_path) {
    int i, j, k, n, dims[ndims];
    float time, *b_field[ndims];
    char full_name[1001], outstr[1001];
    
    // Load the magnetic field data and its dimensions
    sprintf(full_name,"HYM_%0.3d.silo",cycle);
    time = ReadTime_SILO(silo_path,full_name,"HYM_mesh",stopmsg);
    ReadVector_SILO(silo_path,full_name,"b_field",b_field,dims,stopmsg);
    
    // Index and interpolation array initializations:
    int zInd[zLen], rInd[rLen], pInd[pLen];
    float zInterp[zLen], rInterp[rLen], pInterp[pLen];
    
    // Index mapping and printing calls:
    Map_Indices(zLoc,zLen,dims[0]-1,zInd,zInterp);
    Map_Indices(rLoc,rLen,dims[1]-1,rInd,rInterp);
    Map_Indices(pLoc,pLen,dims[2],pInd,pInterp);
    
    // Print out the dimensions and the mapped indices:
    // Print_Dims(dims);
    // Print_Indices(zLen,"z",zInd,zInterp);
    // Print_Indices(rLen,"r",rInd,rInterp);
    // Print_Indices(pLen,"p",pInd,pInterp);

    //------------------------------------------------------------------------//
    // Hack overwrite of zInd array for gathering shifted probe data:
    int ipr_east = 125;
    zInd[0] = ipr_east;
    zInd[1] = 256;
    zInd[2] = 512-ipr_east;
    //------------------------------------------------------------------------//
    
    // Write the magnetic field values at the various probe locations
    ofstream file;
    sprintf(full_name,"Probes_%0.3d.dat",cycle);
    OpenOutputFile(file,probe_path,full_name,stopmsg);
    file << "time= " << time << endl;
    for(i=0; i<zLen; i++) {
        for(j=0; j<rLen; j++) {
            for(k=0; k<pLen; k++) {
                n = fn(zInd[i],rInd[j],pInd[k],dims[0],dims[1]);
                sprintf(outstr,"%6.2f  %5.2f  %8.6f  %13.6E  %13.6E  %13.6E\n",
                        61.0*zLoc[i]-30.5,20.3*rLoc[j],2.0*pi*pLoc[k],
                        b_field[0][n],b_field[1][n],b_field[2][n]);
                file << outstr;
            }
        }
    }
    file.close();
    
    for(int m=0; m<ndims; m++)
        delete [] b_field[m];
    
    cout << "    Output: \"" << probe_path << full_name << "\"\n";
}

//============================================================================//
void Map_Indices(float *loc, int len, int dim, int *ind, float *interp) {
    // Find each index and interpolation mapping:
    for(int i=0; i<len; i++) {
        // Round to the nearest integer:
        ind[i] = (int)(round(loc[i]*dim));
        // Protect against using the outer shell index:
        if(ind[i] == dim)
            ind[i]--;
        // Determine the fractional residual about the nearest integer for
        // interpolation purposes (range = [-0.5,+0.5]):
        interp[i] = loc[i]*dim - ind[i];
    }
}

//============================================================================//
void Print_Indices(int len, char *vname, int *ind, float *interp) {
    // Nicely format the printing of a given index/interpolation array pairing
    char ind_str[1001] = "", interp_str[1001] = "", dum_str1[10], dum_str2[10];
    for(int i=0; i<len; i++) {
        if(strcmp(ind_str,"") == 0) {
            sprintf(dum_str1,"%6d",ind[i]);
            sprintf(dum_str2,"%6.3f",interp[i]);
        }
        else {
            sprintf(dum_str1,",%6d",ind[i]);
            sprintf(dum_str2,",%6.3f",interp[i]);
        }
        strcat(ind_str,dum_str1);
        strcat(interp_str,dum_str2);
    }
    printf("     %sInd = {%s}",vname,ind_str);
    printf("\n  %sInterp = {%s}",vname,interp_str);
    printf("\n\n");
}

//============================================================================//
void Reset_pLoc_Interval() {
    // Function to reset the interval of phi to [0,2*pi]
    float phi_frac;
    for(int k=0; k<pLen; k++) {
        phi_frac = pLoc[k] - floor(pLoc[k]);
        if(phi_frac < 0)
            phi_frac = 1.0 + phi_frac;
        pLoc[k] = phi_frac;
    }
}

//============================================================================//
//============================================================================//
