//============================================================================//
/*

Clayton Myers
ASCII_Write.cpp
Created:  18 August 2009
Modified: 17 September 2009

Support functions for writing data to HYM ASCII files.

*/
//============================================================================//
//============================================================================//

#include <HYM_SILO.hpp>

//============================================================================//
//============================================================================//
void WriteArray_ASCII(ofstream&,int,float*,int);
void MakeCoordArray(float*&,float**,int*);
char *MakeFormatStr(char*,int,int);

//============================================================================//
void WriteScalar_ASCII(char *path, char *fname, char vchar, int *dims,
                       float **mesh_coords, double time, float *var, 
                       char *stopmsg) {
    int Nq = dims[0], Nr = dims[1], Ns = dims[2];
    float *coord_arr;
    char outstr[1001];
    ofstream file;
    
    // Open the output file:    
    OpenOutputFile(file,path,fname,stopmsg);
    
    // Write the header strings:
    sprintf(outstr,"scalar\n%c\nt=%7.1f\n%3d %3d %3d\n",vchar,time,Nq,Nr,Ns);
    file << outstr;
    
    // Write the data the file
    MakeCoordArray(coord_arr,mesh_coords,dims);
    WriteArray_ASCII(file,6,coord_arr,Nq+Nr+Ns);
    WriteArray_ASCII(file,6,var,Nq*Nr*Ns);
    
    file.close();
    delete [] coord_arr;
    
    printf("      ASCII file \"%s\" written to disk.\n",fname);
}

//============================================================================//
void WriteVector_ASCII(char *path, char *fname, char vchar, int *dims, 
                       float **mesh_coords, double time, float **vec, 
                       char *stopmsg) {
    int Nq = dims[0], Nr = dims[1], Ns = dims[2];
    float *coord_arr;
    char outstr[1001];
    ofstream file;
    
    // Open the output file:    
    OpenOutputFile(file,path,fname,stopmsg);
    
    // Write the header strings:
    sprintf(outstr,"vector\n%c\nt=%7.1f\n%3d %3d %3d\n",vchar,time,Nq,Nr,Ns);
    file << outstr;
    
    // Write the data the file
    MakeCoordArray(coord_arr,mesh_coords,dims);
    WriteArray_ASCII(file,6,coord_arr,Nq+Nr+Ns);
    WriteArray_ASCII(file,6,vec[0],Nq*Nr*Ns);
    WriteArray_ASCII(file,6,vec[1],Nq*Nr*Ns);
    WriteArray_ASCII(file,6,vec[2],Nq*Nr*Ns);
        
    file.close();
    delete [] coord_arr;
    
    printf("      ASCII file \"%s\" written to disk.\n",fname);
}

//============================================================================//
void MakeCoordArray(float *&coord_arr, float **mesh_coords, int *dims) {
    int n;
    int Nq = dims[0], Nr = dims[1], Ns = dims[2];
    
    coord_arr = new float[Nq+Nr+(Ns-2*Nghost-1)];    
    n = 0;
    for(int i=0; i<Nq; i++) {
        coord_arr[n] = mesh_coords[0][i];
        n++;
    }
    for(int j=0; j<Nr; j++) {
        coord_arr[n] = mesh_coords[1][j];
        n++;
    }
    for(int k=Nghost; k<(Ns-Nghost-1); k++) {
        coord_arr[n] = mesh_coords[2][k];
        n++;
    }
}

//============================================================================//
void WriteArray_ASCII(ofstream &file, int columns, float *data, int length) {
    int Nrows;
    char outstr[1001], *columnstr, *fmtstr;    
    
    Nrows = length/columns;
    fmtstr = " %13.6E";
    columnstr = MakeFormatStr(fmtstr,7,columns);
    for(int nrow=0; nrow<Nrows; nrow++) {
        int n = columns*nrow;
        sprintf(outstr,columnstr,data[n],data[n+1],data[n+2],data[n+3],
                data[n+4],data[n+5]);
        file << outstr;
    }
    for(int n=0; n<length%columns; n++) {
        sprintf(outstr,fmtstr,data[columns*Nrows+n]);
        file << outstr;
    }
    if(length%columns != 0)
        file << endl;

    delete [] columnstr;
}

//============================================================================//
char* MakeFormatStr(char *fmtstr, int fmtlen, int fmtnumber) {
    char *outstr = new char[fmtlen*fmtnumber+2];
    int nout = 0;
    for(int n=0; n<fmtnumber; n++) {
        for(int m=0; m<fmtlen; m++) {
            outstr[nout] = fmtstr[m];
            nout++;
        }
    }
    outstr[nout] = '\n';
    outstr[nout+1] = '\0';
    return outstr;
}

//============================================================================//
//############################################################################//
// Functional legacy code for reading from single-cycle HYM ASCII files.
//############################################################################//
//============================================================================//
/*
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
*/
//============================================================================//
//============================================================================//
