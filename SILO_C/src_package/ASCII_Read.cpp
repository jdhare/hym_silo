//============================================================================//
/*

Clayton Myers
ASCII_Read.cpp
Created:  17 September 2009
Modified: 17 September 2009

Support functions for reading data from HYM ASCII files.

*/
//============================================================================//
//============================================================================//

#include <HYM_SILO.hpp>

//============================================================================//
//============================================================================//
void ScanMesh_ASCII(ifstream&);
void ReadVar_ASCII(ifstream&,float*&,int*);

//============================================================================//
void ReadMesh_ASCII(char* path, char *fname, int* dims, float **mesh_coords,
                    double &time, char *stopmsg) {
    // This function loads the HYM mesh from a data files
    int Nq, Nr, Ns;
    char str[1024];
    ifstream file;
    
    // Open the designated file
    OpenInputFile(file,path,fname,"ascii",stopmsg);
    // Scroll through the header strings
    for(int i=0; i<3; i++)
        file >> str;
        
    // Get the time and dimension numbers
    file >> time;
    file >> Nq;
    file >> Nr;
    file >> Ns;
    
    // Get the z and r mappings
    float *q = new float[Nq];
    float *r = new float[Nr];
    float *s = new float[Ns];
    
    for(int i=0; i<Nq; i++)
        file >> q[i];
    for(int j=0; j<Nr; j++)
        file >> r[j];
    for(int k=0; k<Ns; k++)
        file >> s[k];    
    
    dims[0] = Nq;
    dims[1] = Nr;
    dims[2] = Ns;
    
    mesh_coords[0] = q;
    mesh_coords[1] = r;
    mesh_coords[2] = s;

    file.close();
}

//============================================================================//
void ReadScalar_ASCII(char *path, char *fname, float *&var, int *dims,
                      char *stopmsg) {
    // This function reads a scalar variable into the array var
    ifstream file;
    OpenInputFile(file,path,fname,"ascii",stopmsg);
    // OpenInputFile(fp,path,fname,"r",stopmsg);
    ScanMesh_ASCII(file);    
    ReadVar_ASCII(file, var, dims);
    file.close();
}

//============================================================================//
void ReadVector_ASCII(char *path, char *fname, float **vec, int *dims, 
                      char *stopmsg) {
    ifstream file;
    OpenInputFile(file,path,fname,"ascii",stopmsg);
    ScanMesh_ASCII(file);        
    ReadVar_ASCII(file,vec[0],dims);
    ReadVar_ASCII(file,vec[1],dims);
    ReadVar_ASCII(file,vec[2],dims);
    file.close();
}

//============================================================================//
void ScanMesh_ASCII(ifstream &file) {
    // This simply scans through the lines containing the HYM mesh in the file
    // given by the file pointer fp.
    int i, j, k, Nq, Nr, Ns;
    char str[1024];

    // Scroll through the header strings
    for(int i=0; i<4; i++)
        file >> str;
    file >> Nq;
    file >> Nr;
    file >> Ns;
    
    // Scroll through the z,r,phi mappings
    for(i=0; i<Nq; i++)
        file >> str;
    for(j=0; j<Nr; j++)
        file >> str;
    for(k=0; k<Ns; k++) 
        file >> str;
}

//============================================================================//
void ReadVar_ASCII(ifstream &file, float *&var, int* dims) {
    // This functions reads a 3D scalar variable from the file given by fp
    int Nq = dims[0], Nr = dims[1], Ns = dims[2];
    var = new float[Nq*Nr*Ns];
    int n = 0;
    for(int k=0; k<Ns; k++) {
        for(int j=0; j<Nr; j++) {
            for(int i=0; i<Nq; i++) {
                file >> var[n];
                n++;
            }
        }
    }
}
//============================================================================//
//============================================================================//
