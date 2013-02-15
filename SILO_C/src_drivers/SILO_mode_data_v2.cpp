//============================================================================//
/*

Clayton Myers
SILO_mode_data_v2.cpp
Created:  15 October 2010
          25 July 2011

Program for reading magnetics data for simulated SSX probe diagnostics from
HYM's .silo storage format.

The goal of this version is to directly output the axial profile of the four  
radially-averaged mode energy densities (n=0/1 and pol/tor) for comparison to 
SSX experimental data.  One file is output for each cycle with the four profiles
 as a function of z.  All of the processing details (radial averaging and RCC 
corrections) are now folded into the C++ code.

The primary motivation for these changes is to enable the user to handle 
arbitrary axial locations for the off-midplane probes.  The physics interest 
here is to pursue the hypothesis that the SSX doublet plasma stretches back into 
the coaxial guns at either end of the device.

----------------------------------------

Command line arguments:
    (1) silo_path (string) -- path to the location of the silo database from 
            which the probe data is being extracted.
    (2) out_path (string) -- path to the location of the output file for the
            extracted probe data.
    (3) cycle (int) -- cycle number for extracted data.  If cycle == 0, data 
            will be extracted from all available cycles.
            
*/
//============================================================================//
//============================================================================//

#include <HYM_SILO.hpp>
#include <SILO_Read.hpp>
#include <Integ_Functions.hpp>

//============================================================================//
//============================================================================//
void ReadArgs(int,char**,char*&,char*&,int&,int&);
void Load_RCC_Data(int*,float**&);
void Write_Mode_Data_Profiles(int,char*,char*,float**);
float Compute_Radial_Average(int,float*);

char *silo_name = "HYM";
char *stopmsg = "Stopping SILO mode data profile extraction.";

//============================================================================//
int main(int argc, char *argv[]) {
    int start_cyc, end_cyc, dims[ndims];
    double prev_time = -1.0, time;
    float **wm0_pol_RCC=NULL;
    char *silo_path=NULL, *out_path=NULL, full_name[1001];
    
    // Read the command line arguments:
    ReadArgs(argc,argv,silo_path,out_path,start_cyc,end_cyc);
    
    // Read and format the RCC data for use in the energy density calculations:
    Get_Mesh_Dims(silo_path,"HYM_001.silo","HYM_mesh",dims,stopmsg);
    Load_RCC_Data(dims,wm0_pol_RCC);

    // Write (to a file) the mode data for the requested cycle argument:
    for(int cyc=start_cyc; cyc<=end_cyc; cyc++) {
        sprintf(full_name,"HYM_%0.3d.silo",cyc);
        time = ReadTime_SILO(silo_path,full_name,"HYM_mesh",stopmsg);
        if(time != prev_time) {
            Write_Mode_Data_Profiles(cyc,silo_path,out_path,wm0_pol_RCC);
            prev_time = time;
        }
    }
    
    // Clean up the RCC data array:
    for(int i=0; i<dims[0]; i++)
        delete [] wm0_pol_RCC[i];
    delete [] wm0_pol_RCC;
}

//============================================================================//
//============================================================================//
void ReadArgs(int argc, char **argv, char *&silo_path, char *&out_path,
              int &start_cyc, int &end_cyc) {
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
    out_path = argv[2];
    VerifyPath(silo_path,stopmsg);
    VerifyPath(out_path,stopmsg);
    ConvertToInt(argv[3],start_cyc,stopmsg);
    ConvertToInt(argv[4],end_cyc,stopmsg);
}

//============================================================================//
void Load_RCC_Data(int *dims, float **&wm0_pol_RCC) {
    // Reads the RCC data produced by the vacuum simulation
    // Hard coded for a 257x129 grid from the vacuum simulation
    // Handles either 257x129 or 513x129 from the plasma simulation
    // Error checking for the above criterion is not robust.
    
    char *RCC_path = "./RCC_Data/"; 
    char *RCC_fname = "SILO_mode_data_RCC_257_129_cyc265.dat";

    wm0_pol_RCC = new float*[dims[0]];

    // Read the RCC data to the local storage array:
    ifstream file;
    OpenInputFile(file,RCC_path,RCC_fname,"ascii",stopmsg);
    for(int i=0; i<dims[0]; i++) {
        wm0_pol_RCC[i] = new float[dims[1]];
        if(dims[0]==257 || (dims[0]==513 && i%2==0)) {
            float dum;  file >> dum;  file >> dum;
            for(int j=0; j<dims[1]; j++)
                file >> wm0_pol_RCC[i][j];
        }
    }
    file.close();

    // Interpolate from Nz=257 to Nz=513 if necessary:
    if(dims[0]==513) {
        for(int i=1; i<dims[0]; i+=2) {
            for(int j=0; j<dims[1]; j++)
                wm0_pol_RCC[i][j] = (wm0_pol_RCC[i+1][j]+wm0_pol_RCC[i-1][j])/2.;
        }
    }
}

//============================================================================//
void Write_Mode_Data_Profiles(int cycle, char *silo_path, char *out_path,
                              float **wm0_pol_RCC) {
    int i, j, n, m, dims[ndims];
    float time, *b_field[ndims];
    char full_name[1001], outstr[1001];
    
    //------------------------------------------------------------------------//
    // Load the magnetic field data and its dimensions:
    sprintf(full_name,"HYM_%0.3d.silo",cycle);
    time = ReadTime_SILO(silo_path,full_name,"HYM_mesh",stopmsg);
    ReadVector_SILO(silo_path,full_name,"b_field",b_field,dims,stopmsg);
    
    //------------------------------------------------------------------------//
    // Initialize the Fourier coefficients c0 and c1:
    float *c0[ndims], *c1[ndims];
    for(m=0; m<ndims; m++) {
        c0[m] = new float[dims[0]*dims[1]*dims[2]];
        c1[m] = new float[dims[0]*dims[1]*dims[2]];
    }
    
    //------------------------------------------------------------------------//
    // Compute the full z-r array of Fourier-decomposed Brms coefficients:
    for(i=0; i<dims[0]; i++) {
        for(j=0; j<dims[1]; j++) {
            n = fn(i,j,0,dims[0],dims[1]);
            for(m=0; m<ndims; m++) {
                Fourier_Decomp(b_field,dims,i,j,m,c0[m][n],c1[m][n]);
            }
        }
    }

    //------------------------------------------------------------------------//
    // Compute and save the mode energy densities:
    int Nr = dims[1];
    float wm0_pol_rad[Nr], wm0_tor_rad[Nr], wm1_pol_rad[Nr], wm1_tor_rad[Nr];
    float wm0_pol_avg, wm0_tor_avg, wm1_pol_avg, wm1_tor_avg;

    // Open this cycle's data file for saving the energy density profiles:
    ofstream file;
    sprintf(full_name,"SILO_mode_data_axial_profiles_%03d.dat",cycle);
    OpenOutputFile(file,out_path,full_name,stopmsg);
    sprintf(outstr,"%6.1f",time);
    file << "time= " << outstr << endl;
    
    for(i=0; i<dims[0]; i++) {
        // First compute the radial profiles at this z location:
        for (j=0; j<dims[1]; j++) {
            n = fn(i,j,0,dims[0],dims[1]);
            wm0_pol_rad[j] = c0[0][n]*c0[0][n] + c0[1][n]*c0[1][n];
            wm0_tor_rad[j] = c0[2][n]*c0[2][n];
            wm1_pol_rad[j] = c1[0][n]*c1[0][n] + c1[1][n]*c1[1][n];
            wm1_tor_rad[j] = c1[2][n]*c1[2][n];

            // Subtract the RCC energy from the wm0_pol radial profile:
            wm0_pol_rad[j] -= wm0_pol_RCC[i][j];
        }

        // Compute the radial averages:
        wm0_pol_avg = Compute_Radial_Average(Nr,wm0_pol_rad);
        wm0_tor_avg = Compute_Radial_Average(Nr,wm0_tor_rad);
        wm1_pol_avg = Compute_Radial_Average(Nr,wm1_pol_rad);
        wm1_tor_avg = Compute_Radial_Average(Nr,wm1_tor_rad);

        // Write the results to the data file:
        sprintf(outstr,"%03d%8.2f",i,i*61.0/(1.*dims[0]-1)-30.5);
        file << outstr;
        sprintf(outstr,"%14.6E%14.6E",wm0_pol_avg,wm0_tor_avg);
        file << outstr;
        sprintf(outstr,"%14.6E%14.6E",wm1_pol_avg,wm1_tor_avg);
        file << outstr << endl;
    }

    file.close();
    
    //------------------------------------------------------------------------//
    for(m=0; m<ndims; m++) {
        delete [] b_field[m];
        delete [] c0[m];
        delete [] c1[m];
    }        

    printf("    Completed Cycle %03d\n",cycle);
}

//============================================================================//
float Compute_Radial_Average(int Nr, float *wmn_rad) {
    float rj_mid, wmn_mid, sum_wmn=0.0;

    float Rc = 20.3;
    float dr = Rc/(Nr-1.);

    for (int j=1; j<Nr; j++) {
        rj_mid = dr * (j-1/2.);
        wmn_mid = (wmn_rad[j]+wmn_rad[j-1])/2.;
        sum_wmn += 2.*pi*dr*rj_mid*wmn_mid;
    }
    return sum_wmn/(pi*Rc*Rc);
}

//============================================================================//
//============================================================================//
