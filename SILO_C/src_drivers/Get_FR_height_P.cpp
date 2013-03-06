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

char *stopmsg   = "Stopping FR height extraction.";

//============================================================================//
int main(int argc, char *argv[]) {

    char *silo_path=NULL, *out_path=NULL, fname[1001];
    
    // Read the command line arguments:
    ReadArgs(argc,argv,silo_path,out_path);
    
    //get number of silo cycles
    int Ncyc = Get_Ncyc(silo_path, stopmsg);
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
}

//============================================================================//
void Get_Max_Vals(int cycle, char *silo_path, FILE *fp) {
    int dims[ndims], ib1, ib2, jb1, jb2, i, j, k, n, z_loc;
    float *pressure;
    float p;
    double time;
    char full_name[1001];
    
    // Load the magnetic field data and its dimensions:
    sprintf(full_name,"HYM_%0.3d.silo",cycle);
    time = ReadTime_SILO(silo_path,full_name,"HYM_mesh",stopmsg);
    ReadScalar_SILO(silo_path,full_name,"pressure",pressure,dims,stopmsg);

    
    // Set the boundary indices:
    ib1 = 0; //   3;
    ib2 = dims[0]; // 509;
    jb1 =  0; //   1;
    jb2 = dims[1]; // 118;
    k=dims[2]/2; //Look only at midplane
    
    // Initialize the maxJ and maxv variables:
    float p_max=0.0;
    
    // Gather Jmax and vmax for this cycle:
    for(i=ib1; i<ib2; i++) {
        for(j=jb1; j<jb2; j++) {
            p=0.0; z_loc=0.0;
            n = fn(i,j,k,dims[0],dims[1]);
            p=pressure[n];
            if(p > p_max){
                p_max=p;
                z_loc=j;
            }
        }
    }
    
    delete pressure;

    //------------------------------------------------------------------------//
    // Write out the max values for the cycle:

    fprintf(fp,"%03d%12.4f%12.4f%03d\n",
               cycle,time,p_max,z_loc);
 
    printf("    Completed Cycle %03d\n",cycle);
}

//============================================================================//
//============================================================================//

