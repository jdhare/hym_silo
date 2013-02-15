//============================================================================//
/*

Clayton Myers
SILO_Read.cpp
Created:  14 September 2009
Modified: 16 September 2009

Support functions for reading data from existing .silo databases.

*/
//============================================================================//
//============================================================================//

#include <HYM_SILO.hpp>
#include <SILO_Read.hpp>

//============================================================================//
//============================================================================//
void GetVar_SILO(DBfile*,char*,DBquadvar*&,char*,int*,char*);
void StripVar_SILO(DBquadvar*,float*&,int*,int);
void Cart_to_Cyl(float**,int*);

//============================================================================//
void ReadScalar_SILO(char *path, char *fname, char *varname, float *&var, 
                     int *dims, char *stopmsg) {
    // Function to read a scalar variable from an existing .silo database
    DBfile *dbfile=NULL;
    OpenFile_SILO(path,fname,dbfile,stopmsg);
    
    DBquadvar *dbvar=NULL;
    GetVar_SILO(dbfile,varname,dbvar,"scalar",dims,stopmsg);        
    StripVar_SILO(dbvar,var,dims,0);
     
    DBFreeQuadvar(dbvar);
    DBClose(dbfile);
}

//============================================================================//
void ReadVector_SILO(char *path, char *fname, char *varname, float **vec,
                     int *dims, char *stopmsg) {
    // Function to read a vector variable from an existing .silo database
    DBfile *dbfile=NULL;
    OpenFile_SILO(path,fname,dbfile,stopmsg);
    
    DBquadvar *dbvar=NULL;
    GetVar_SILO(dbfile,varname,dbvar,"vector",dims,stopmsg);        
    StripVar_SILO(dbvar,vec[0],dims,0);
    StripVar_SILO(dbvar,vec[1],dims,1);
    StripVar_SILO(dbvar,vec[2],dims,2);
    
    Cart_to_Cyl(vec,dims);
    
    DBFreeQuadvar(dbvar);
    DBClose(dbfile);                     
}

//============================================================================//
double ReadTime_SILO(char *path, char *fname, char *mesh_name, char *stopmsg) {
    DBfile *dbfile=NULL;
    OpenFile_SILO(path,fname,dbfile,stopmsg);
    double time;
        
    DBquadmesh *dbmesh=NULL;
    dbmesh = DBGetQuadmesh(dbfile,mesh_name);
    
    time = dbmesh->dtime;
    
    DBFreeQuadmesh(dbmesh);
    DBClose(dbfile);
    
    return time;
}

//============================================================================//
//============================================================================//
void GetVar_SILO(DBfile *dbfile, char *varname, DBquadvar *&dbvar, 
                 char *nvals_str, int *dims, char *stopmsg) {
    // Function to retrieve variable from .silo database with error checking
    char message[1001];    
    // Check whether the variable requested exists in the .silo database
    DBtoc *dbtoc = DBGetToc(dbfile);
    bool present = false;
    for(int i=0; i<dbtoc->nqvar; i++) {
        if(strcmp(varname,dbtoc->qvar_names[i]) == 0)
            present = true;
    }
    // Stop execution if the variable requested does not exist in the database
    if(!present) {
        sprintf(message,"  The variable \"%s\" %s\n  %s",varname,
                "was not found in the .silo database.",stopmsg);
        StopExecution(message);
    }
    // Get the variable
    dbvar = DBGetQuadvar(dbfile,varname);    
    // Stop execution if a variable of improper rank is requested
    int nvals;
    if(strcmp(nvals_str,"scalar") == 0)
        nvals = 1;
    else if(strcmp(nvals_str,"vector") == 0)
        nvals = 3;
    else
        nvals = 0;
    if(dbvar->nvals != nvals) {
        sprintf(message,"  The variable \"%s\" %s\n  %s%d\n  %s%d\n  %s",
                varname,"does not have the proper rank.","  Requested rank: ",
                nvals,"  Variable rank:  ",dbvar->nvals,stopmsg);
        StopExecution(message);    
    } 
    // Stop execution if the variable access call fails
    if(dbvar == NULL) {
        sprintf(message,"  %s \"%s\" %s\n  %s",
                "The call to access the variable",varname,
                "failed for an unknown reason.",stopmsg);
        StopExecution(message);
    }
    
    // Transfer the dimension array values
    for(int m=0; m<ndims; m++)
        dims[m] = dbvar->dims[m];
}

//============================================================================//
void StripVar_SILO(DBquadvar *dbvar, float *&var, int *dims, int component) {
    // Function to strip the ghost zones from the SILO variable
    int i,j,k,n1,n2,Ntot,min_phi,max_phi;

    // Calculate the stripped array dimension and allocate the new array
    min_phi = dbvar->min_index[2]; // Lower ghost zone index
    max_phi = dbvar->max_index[2]; // Upper ghost zone index
    dims[2] = max_phi-min_phi-1;
    Ntot = dims[0]*dims[1]*dims[2];
    var = new float[Ntot];
    
    // Transfer the values between the arrays
    // n1: counter for stripped array
    // n2: counter for unstripped array
    n1 = 0;
    for(k=min_phi; k<(max_phi-1); k++) {
        for(j=0; j<dims[1]; j++) {
            for(i=0; i<dims[0]; i++) {
                n2 = fn(i,j,k,dims[0],dims[1]);
                var[n1] = dbvar->vals[component][n2];
                n1++;
            }
        }
    }
}

//============================================================================//
void Cart_to_Cyl(float **vec, int *dims) {
    // Add functionality comment here
    int i, j, k, n;
    int Nq = dims[0], Nr = dims[1], Ns = dims[2];
    int Ntot = Nq*Nr*Ns;
    float *vec_x = vec[0], *vec_y = vec[1], *vec_q = vec[2], *s=NULL;
    float *vec_r = new float[Ntot];
    float *vec_s = new float[Ntot];
    
    // Construct the phi vector artificially:
    Construct_Phi(s,dims);
    
    // Convert (x,y) components to (r,s) components
    n = 0;
    for(k=0; k<Ns; k++) {
        for(j=0; j<Nr; j++) {
            for(i=0; i<Nq; i++) {
                if(j == 0) {
                    vec_r[n] = 0.0;
                    vec_s[n] = 0.0;
                }
                else {
                    vec_r[n] =  vec_x[n]*cos(s[k]) + vec_y[n]*sin(s[k]);
                    vec_s[n] = -vec_x[n]*sin(s[k]) + vec_y[n]*cos(s[k]);
                }
                n++;
            }
        }
    }
    
    // Replace and sort the vector components into the (q,r,s) arrangement
    delete [] vec_x; delete [] vec_y; delete [] s;
    vec[0] = vec_q;
    vec[1] = vec_r;
    vec[2] = vec_s;
}

//============================================================================//
//============================================================================//
void OpenFile_SILO(char *path, char *fname, DBfile *&dbfile, char *stopmsg) {
    // Open an existing SILO file to read data
    char fullpath[1001];
    strcpy(fullpath,path);
    strcat(fullpath,fname);
    dbfile = NULL;
    dbfile = DBOpen(fullpath,DB_PDB,DB_APPEND);
    if(dbfile == NULL) {
        char message[1001];
        sprintf(message,"  The SILO database \"%s\" %s\n  %s",fullpath,
                "could not be opened.",stopmsg);
        StopExecution(message);
    }
}

//============================================================================//
int Get_Ncyc(char *silo_path, char *stopmsg) {
    // Function to find the number of available cycles based on the files in the
    // silo path directory.
    int nlist, Ncyc, Ncyc_old;
    char **silo_list, testname[1001];
    
    // Get a directory listing of the silo path:
    listdir(silo_path,silo_list,nlist,stopmsg);
    
    // Iterate through the silo path to find which files exist:
    Ncyc = 0;
    do {
        Ncyc_old = Ncyc;
        sprintf(testname,"HYM_%03d.silo",Ncyc+1);
        for(int m=0; m<nlist; m++) {
            // If the file is in the list, increment Ncyc:
            if(strcmp(silo_list[m],testname) == 0)
                Ncyc++;
        }
    } while(Ncyc_old != Ncyc); // Go until the file testname doesn't exist

    delete [] silo_list;
    
    if(Ncyc == 0) {
        char message[1001];
        sprintf(message,"  %s%s%s\n  %s",
                "The directory \"",silo_path,
                "\" contains no source SILO databases.",stopmsg);
        StopExecution(message);
    }
    
    return Ncyc;
}

//============================================================================//
void Get_Mesh_Dims(char *path, char *fname, char *mesh_name, int *dims, 
                   char *stopmsg) {
    // Function to read the dimensions of the mesh in an existing .silo database
    char message[1001];   
    DBfile *dbfile=NULL;
    OpenFile_SILO(path,fname,dbfile,stopmsg);
    
    DBquadmesh *dbquadmesh=NULL;
    dbquadmesh = DBGetQuadmesh(dbfile,mesh_name);
    if(dbquadmesh == NULL) {
        sprintf(message,"  Unable to locate the mesh \"%s\" %s\n  %s",mesh_name,
                "in the .silo database.",stopmsg);
        StopExecution(message);
    }
    
    for(int m=0; m<ndims; m++)
        dims[m] = dbquadmesh->dims[m];
    dims[2] = dims[2] - (2*Nghost+1);
    
    DBFreeQuadmesh(dbquadmesh);
    DBClose(dbfile);
}

//============================================================================//
//============================================================================//
