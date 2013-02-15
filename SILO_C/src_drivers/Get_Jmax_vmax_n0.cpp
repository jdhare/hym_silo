//============================================================================//
/*

Clayton Myers
Lambda_Calc.cpp
Created:  05 February 2010

Command line arguments:
    (1) silo_path (string) -- path to the location of the silo database(s).
    (2) out_path (string) -- path to the location of the output file for the
            calculated lambda data.
    (3) cycle (int) -- cycle number for extracted data.  If cycle == 0, data 
            will be extracted from all available cycles.

*/
//============================================================================//
//============================================================================//

#include <HYM_SILO.hpp>
#include <SILO_Read.hpp>

#define WRITE false

//============================================================================//
//============================================================================//
void ReadArgs(int,char**,char*&,char*&);
void Get_Max_Vals(int,char*,FILE*);

char *stopmsg   = "Stopping Jmax/vmax extraction.";

//============================================================================//
int main(int argc, char *argv[]) {
    int Ncyc = 210;
    char *silo_path=NULL, *out_path=NULL, fname[1001];
    
    // Read the command line arguments:
    ReadArgs(argc,argv,silo_path,out_path);
    
    // Open the file for saving the data:
    sprintf(fname,"%s/Jmax_vmax_n0_w_mins.dat",out_path);
    FILE *fp = fopen(fname,"w");
    
    // Write (to a file) the probe data for the requested cycle argument:
    for(int cyc=1; cyc<=Ncyc; cyc++)
        Get_Max_Vals(cyc,silo_path,fp);
        
    fclose(fp);
}

//============================================================================//
//============================================================================//
void ReadArgs(int argc, char **argv, char *&silo_path, char *&out_path) {
    // Count the initial command line arguments:
    if(!(argc == (1+2))) {
        char message[1001];
        sprintf(message,"  %s\n  %s",
                "An improper number of command line arguments was found.",
                stopmsg);
        StopExecution(message);
    }
        
    // Distribute the command line arguments:
    silo_path = argv[1];
    out_path = argv[2];
    VerifyPath(silo_path,stopmsg);
    VerifyPath(out_path,stopmsg);
    // ConvertToInt(argv[3],cycle,stopmsg);
}

//============================================================================//
void Get_Max_Vals(int cycle, char *silo_path, FILE *fp) {
    int dims[ndims], ib1, ib2, jb1, jb2, i, j, k, n;
    float *current[ndims], *velocity[ndims];
    float J[3], v[3], Jpol, Jtor, vpol, vtor, vr, vtheta; 
    double time;
    char full_name[1001];
    
    // Load the magnetic field data and its dimensions:
    sprintf(full_name,"HYM_%0.3d.silo",cycle);
    time = ReadTime_SILO(silo_path,full_name,"HYM_mesh",stopmsg);
    ReadVector_SILO(silo_path,full_name,"current_density",current,dims,stopmsg);
    ReadVector_SILO(silo_path,full_name,"velocity",velocity,dims,stopmsg);
    
    // Set the boundary indices:
    ib1 = 200; //   3;
    ib2 = 312; // 509;
    jb1 =  32; //   1;
    jb2 = 120; // 118;
    
    // Initialize the maxJ and maxv variables:
    float Jpol_max = 0.0, Jtor_max = 0.0;
    float vpol_max = 0.0, vtor_max = 0.0;
    float vr_max = 0.0, vr_min = 0.0;
    float vtheta_max = 0.0, vtheta_min = 0.0; 
    
    // Gather Jmax and vmax for this cycle:
    for(i=ib1; i<ib2; i++) {
        for(j=jb1; j<jb2; j++) {
            Jpol = 0.0; Jtor = 0.0;
            vpol = 0.0; vtor = 0.0;        
            for(k=0; k<dims[2]; k++) {
                n = fn(i,j,k,dims[0],dims[1]);
                for(int m=0; m<ndims; m++) {
                    J[m] = current[m][n];
                    v[m] = velocity[m][n];
                }
                Jpol += sqrt(J[0]*J[0] + J[1]*J[1]);
                Jtor += sqrt(J[2]*J[2]);
                vpol += sqrt(v[0]*v[0] + v[1]*v[1]);
                vtor += sqrt(v[2]*v[2]);
                vr += v[1];
                vtheta += v[2];           
            }
            Jpol = Jpol/dims[2]; Jtor = Jtor/dims[2];
            vpol = vpol/dims[2]; vtor = vtor/dims[2];
            vr = vr/dims[2]; vtheta = vtheta/dims[2];

            if(Jpol > Jpol_max)
                Jpol_max = Jpol;
            if(Jtor > Jtor_max)
                Jtor_max = Jtor;
            if(vpol > vpol_max)
                vpol_max = vpol;
            if(vtor > vtor_max)
                vtor_max = vtor;
                
            if(vr < vr_min)
                vr_min = vr;
            if(vr > vr_max)
                vr_max = vr;
            if(vtheta < vtheta_min)
                vtheta_min = vtheta;
            if(vtheta > vtheta_max)
                vtheta_max = vtheta;                
        }
    }
    
    for(int m=0; m<ndims; m++) {
        delete [] current[m];
        delete [] velocity[m];
    }

    //------------------------------------------------------------------------//
    // Write out the max values for the cycle:

    fprintf(fp,"%03d%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f\n",
               cycle,time,Jpol_max,Jtor_max,vpol_max,vtor_max,
               vr_min,vr_max,vtheta_min,vtheta_max);
 
    printf("    Completed Cycle %03d\n",cycle);
}

//============================================================================//
//============================================================================//

