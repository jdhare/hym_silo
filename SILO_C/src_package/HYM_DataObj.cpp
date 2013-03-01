//============================================================================//
/*

Clayton Myers
HYM_DataObj.cpp
Created:  15 September 2009
Modified: 06 February 2010

Class for HYM scalar and vector data objects.

*/
//============================================================================//
//============================================================================//

#include <HYM_SILO.hpp>
#include <HYM_DataObj.hpp>
#include <ASCII_Write.hpp>

//============================================================================//
const int intsize = sizeof(int);
const int dblsize = sizeof(double);

// There are 4 ghost zones (2 on either end) in z (q) and phi (s).  The r = 0
// and r = dr (j = 0 and j = 1) points are not considered ghost zones outside of
// HYM, so there are only 2 stripped ghost zones in r:
const int Nghost_q1 = 2;
const int Nghost_q2 = 2;
const int Nghost_r1 = 0;
const int Nghost_r2 = 2;
const int Nghost_s1 = 2;
const int Nghost_s2 = 2;

//============================================================================//
//############################################################################//
//============================================================================//
HYMDataObj::HYMDataObj(char vchar, char *data_path, char *fname_src,
                       int *dims, char *stopmsg, int nvals) {
    // Assign basic members:
    this->data_path = data_path;
    this->fname_src = fname_src;
    this->vchar = vchar;
    this->stopmsg = stopmsg;
    this->cycle_mask=NULL;
    this->times = NULL;
    
    // Verify the source data file:
    if(!VerifyInputFile(data_path,fname_src,"binary",stopmsg)) {
        char message[1001];
        sprintf(message,"      %s%c%s\n      %s%s\n      %s%s\n      %s",
                "The data for the variable ",vchar," was not found:",
                "    Data Path:    ",data_path,"    Missing File: ",fname_src,
                stopmsg);
        StopExecution(message);
    }
    
    // Open the source data file:
    OpenInputFile(this->file,data_path,this->fname_src,"binary",this->stopmsg);
    
    // Assign mesh and dimension attribute members:
    this->dims = dims;
    this->Ntot = dims[0]*dims[1]*dims[2];
    
    this->dims_in[0] = Nghost_q1 + dims[0] + Nghost_q2;
    this->dims_in[1] = Nghost_r1 + dims[1] + Nghost_r2;
    this->dims_in[2] = Nghost_s1 + dims[2] + Nghost_s2;
    this->Ntot_in = dims_in[0]*dims_in[1]*dims_in[2];
    
    // Assign data attribute members:
    this->nvals = nvals;
    if(nvals == 1)
        this->data_type = "scalar";
    else if(nvals == 3)
        this->data_type = "vector";
    else
        StopExecution("  Invalid nvals argument in HYMDataObj constructor.");
    this->record_length = (long)(5*intsize + (1+Ntot_in*nvals)*dblsize + 20);
    
    // Validate the source data file and build the time vector:
    this->ValidateSourceFile();
    this->GetTimes();
}

//============================================================================//
HYMDataObj::~HYMDataObj(void) { 
    this->file.close();
    if(this->cycle_mask != NULL)
        delete [] this->cycle_mask;
    if(this->times != NULL)
        delete [] this->times;
}

//============================================================================//
//============================================================================//
void HYMDataObj::ValidateSourceFile(void) {
    // This function examines a source file to determine how many complete data
    // cycles are stored in the file.
    int dims1[ndims], dims2[ndims], Ncyc;
    long file_length, cycle_pos;

    // First determine the number of cycles that are stored from the file size:
    this->file.seekg(0,ios::end);
    file_length = this->file.tellg();
    Ncyc = 0;
    while((Ncyc*this->record_length) < file_length)
        Ncyc++;

    // Determine which cycles have stored data and save the result to the
    // cycle_mask data member:        
    if(Ncyc > 0)
        this->cycle_mask = new bool[Ncyc];
    else {
        char message[1001];
        sprintf(message,"  %s%s%s\n  %s\n      %s",
                "    Error validating the binary file \"",this->fname_src,"\":",
                "        No data is stored in the file.",stopmsg);
        StopExecution(message);
    }
    this->Ncyc = Ncyc;
    
    // Validate the stored data for each cycle:
    for(int m=0; m<Ncyc; m++) {
        bool valid = false;
        // First make sure the entire cycle is within the file length:
        cycle_pos = m*this->record_length;
        valid = cycle_pos < file_length;
        // Compare the mesh dimensions of this cycle with the global mesh
        // dimensions:
        if(valid) {
            this->ReadDims(cycle_pos+2*intsize+dblsize,dims1);
            dims1[0] = dims1[0] - Nghost_q2 - Nghost_q1;
            dims1[1] = dims1[1] - Nghost_r2 - Nghost_r1;
            dims1[2] = dims1[2] - Nghost_s2 - Nghost_s1;
            valid = HYMDataObj::CompareMeshDims(this->dims,dims1,
                                                this->fname_src,this->stopmsg);
        }            
        // Check if the entire data record exists for this cycle:
        if(valid)
            valid = !((cycle_pos+this->record_length) > file_length);        
        // Assign the cycle_mask value:
        this->cycle_mask[m] = valid;
    }

    //------------------------------------------------------------------------//
    // Write out the validation results:
    char fname_str[1001];
    sprintf(fname_str,"%-10.10s",this->fname_src);
    cout << "      Reading: " << fname_str;
    cout << " Ncyc = " << Ncyc << endl;
}

//============================================================================//
void HYMDataObj::ReadDims(long pos, int *dims) {
    this->file.seekg(pos,ios::beg);
    this->file.read((char*)&dims[0],intsize);
    this->file.read((char*)&dims[1],intsize);
    this->file.read((char*)&dims[2],intsize);
}

//============================================================================//
void HYMDataObj::GetTimes(void) {
    // Reads the time at each cycle into a single vector
    this->times = new double[this->Ncyc];    
    for(int cycle=1; cycle<=Ncyc; cycle++) { 
        this->file.seekg((cycle-1)*this->record_length + 2*intsize,ios::beg);
        this->file.read((char*)&this->times[cycle-1],dblsize);
    }
}

//============================================================================//
//============================================================================//
void HYMDataObj::ReadMesh_Binary(char *data_path, char* fname_mesh, int *dims,
                                 float **mesh_coords, char *stopmsg) {
    // This function loads the HYM mesh from a data file (3D or 2D)
    int k, kshift, Nq, Nr, Ns, Ntot;
    char str[1024];
    ifstream file;
    
    // Open the designated file:
    OpenInputFile(file,data_path,fname_mesh,"binary",stopmsg);

    // Scroll through unneeded variables:
    file.seekg(6*intsize,ios::beg);

    // Get the input dimensions:
    file.read((char*)&Nq,intsize);
    file.read((char*)&Nr,intsize);
    file.read((char*)&Ns,intsize);
    
    // Read the mesh coordinates:
    file.seekg(dblsize,ios::cur);
    double *q_buffer = new double[Nq];
    double *r_buffer = new double[Nr];
    double *s_buffer = new double[Ns];
    file.read((char*)q_buffer,Nq*dblsize);
    file.read((char*)r_buffer,Nr*dblsize);
    file.read((char*)s_buffer,Ns*dblsize);
    
    file.close();
    
    // Strip off the HYM ghost coordinates:
    StripGhostCoords(q_buffer,mesh_coords[0],Nq,Nghost_q1,Nghost_q2);
    StripGhostCoords(r_buffer,mesh_coords[1],Nr,Nghost_r1,Nghost_r2);
    StripGhostCoords(s_buffer,mesh_coords[2],Ns,Nghost_s1,Nghost_s2);
    
    // Validate the dims values:
    int fdims[ndims];
    fdims[0] = Nq - Nghost_q2 - Nghost_q1;
    fdims[1] = Nr - Nghost_r2 - Nghost_r1;
    fdims[2] = Ns - Nghost_s2 - Nghost_s1;
    HYMDataObj::CompareMeshDims(fdims,dims,fname_mesh,stopmsg);
}

//============================================================================//
void HYMDataObj::StripGhostCoords(double *buffer, float *&stripped_coords, 
                                  int N, int Nghost_1, int Nghost_2) {
    // Strips the HYM ghost zones off of a single coordinate vector
    stripped_coords = new float[N-Nghost_2-Nghost_1];
    int n = 0;    
    for(int m=Nghost_1; m<(N-Nghost_2); m++) {
        stripped_coords[n] = (float)buffer[m];
        n++;
    }
    delete [] buffer;
}

//============================================================================//
bool HYMDataObj::CompareMeshDims(int *dims1, int *dims2, char *fname, 
                                 char *stopmsg) {
    // Checks if the file mesh has the same dimensions as the source mesh:
    if(dims1[0] != dims2[0] || dims1[1] != dims2[1] || dims1[2] != dims2[2]) {
        return false;
        /*
        char dimstr1[1001], dimstr2[1001], message[1001];
        sprintf(dimstr1,"(%3d,%3d,%3d)",dims1[0],dims1[1],dims1[2]);
        sprintf(dimstr2,"(%3d,%3d,%3d)",dims2[0],dims2[1],dims2[2]);
        sprintf(message,"  %s%s%s\n  %s\n          %s <--> %s\n  %s\n      %s",
                "    Error validating the binary file \"",fname,"\":",
                "    The file does not have the proper mesh dimensions:",
                dimstr1,dimstr2,
                "    Check to make sure the correct file is being accessed.",
                stopmsg);
        StopExecution(message);
        */  
    }
    else
        return true;
}

//============================================================================//
//============================================================================//
void HYMDataObj::PositionPointer_Binary(int cycle) {
    // Error check for the requested cycle number:
    if(cycle < 1 || cycle > this->Ncyc) {
        char message[1001];
        sprintf(message,"      %s%s%s\n      %s%d%s\n      %s%d%s\n      %s",
                "Error accessing the binary file \"",this->fname_src,"\":",
                "The requested cycle number is not valid (cycle = ",cycle,").",
                "The cycle number must be in the range [1,",this->Ncyc,"].",
                stopmsg);
        StopExecution(message);
    }
    // Move the pointer to the head of the data:
    this->file.seekg((cycle-1)*record_length + 5*intsize + dblsize,ios::beg);
}

//============================================================================//
void HYMDataObj::ReadVar_Binary(float *&var) {
    // Reads a single variable from the source file (pointer must be 
    // prepositioned by a function such as PositionPointer_Binary).
    double *var_buffer = new double[this->Ntot_in];
    this->file.read((char*)var_buffer,this->Ntot_in*dblsize);    
        
    // Strip the HYM ghost zones (in z, r, and phi) from the variable.
    var = new float[this->Ntot];
    int n = 0;
    for(int k=Nghost_s1; k<(dims_in[2]-Nghost_s2); k++) {
        for(int j=Nghost_r1; j<(dims_in[1]-Nghost_r2); j++) {
            for(int i=Nghost_q1; i<(dims_in[0]-Nghost_q2); i++) {
                var[n] = (float)var_buffer[fn(i,j,k,dims_in[0],dims_in[1])];
                n++;
            }
        }
    }
    delete [] var_buffer;
}

//============================================================================//
//############################################################################//
//============================================================================//
HYMScalarObj::HYMScalarObj(char vchar, char* data_path, char * fname_src,
                           int *dims, char* stopmsg) :
              HYMDataObj(vchar,data_path,fname_src,dims,stopmsg,1) 
{
    this->id='S';//scalar
    if(vchar == 'p') {
        this->varname = "pressure";
        this->ascii_name = "p3out";
    }
    else if(vchar == 'n') {
        this->varname = "density";
        this->ascii_name = "n3out";
    }
    else {
        char message[1001];
        sprintf(message,"  %s%c%s\n      %s",
                "Unrecognized vchar argument \'",this->vchar,
                "\' in HYMVectorObj constructor",this->stopmsg);
        StopExecution(message);
    }
}

//============================================================================//
 void HYMScalarObj::GetData_SILO(int cycle, float *&var){
     this->ReadScalar_Binary(cycle, var);
 }
//These functions used to be used to write silo databases from inside the HYM data object.
//This functionality (almost exactly) can now be found inside Silo Cyc Obj.
////============================================================================//
//void HYMScalarObj::WriteData_SILO(DBfile *dbfile, int cycle, char *mesh_name,
//                                  float **mesh_coords) {
//    float *var;
//    this->ReadScalar_Binary(cycle,var);
//    WriteScalar_SILO(dbfile,this->varname,mesh_name,var,this->dims);
//    delete [] var;
//}

//============================================================================//
void HYMScalarObj::WriteData_ASCII(char *ascii_path, int cycle, double tout, 
                                   float **mesh_coords) {
    char fname_ascii[1001];
    sprintf(fname_ascii,"%s_%03d.dat",this->ascii_name,cycle);
    float *var;
    this->ReadScalar_Binary(cycle,var);
    WriteScalar_ASCII(ascii_path,fname_ascii,this->vchar,this->dims,   
                      mesh_coords,tout*cycle,var,this->stopmsg);
    delete [] var;
}

//============================================================================//
void HYMScalarObj::ReadScalar_Binary(int cycle, float *&var) {
    this->PositionPointer_Binary(cycle);
    this->ReadVar_Binary(var);
}

//============================================================================//
//############################################################################//
//============================================================================//
HYMVectorObj::HYMVectorObj(char vchar, char* data_path, char * fname_src,
                           int *dims, char* stopmsg) :
              HYMDataObj(vchar,data_path,fname_src,dims,stopmsg,3) 
{
    this->id='V';//vector
    if(vchar == 'B') {
        this->varname = "b_field";
        this->varnames[0] = "B_x"; 
        this->varnames[1] = "B_y";
        this->varnames[2] = "B_z";
        this->ascii_name = "b3out";
    }
    else if(vchar == 'v') {
        this->varname = "velocity";        
        this->varnames[0] = "v_x"; 
        this->varnames[1] = "v_y";
        this->varnames[2] = "v_z";
        this->ascii_name = "v3out";
    }
    else if(vchar == 'J') {
        this->varname = "current_density";
        this->varnames[0] = "J_x"; 
        this->varnames[1] = "J_y";
        this->varnames[2] = "J_z";
        this->ascii_name = "j3out";
    }
    else {
        char message[1001];
        sprintf(message,"  %s%c%s\n      %s",
                "Unrecognized vchar argument \'",this->vchar,
                "\' in HYMVectorObj constructor",this->stopmsg);
        StopExecution(message);
    }
}

//============================================================================//
 void HYMVectorObj::GetData_SILO(int cycle, float **vecs){
     this->ReadVector_Binary(cycle, vecs);
 }
//============================================================================//
//These functions used to be used to write silo databases from inside the HYM data object.
//This functionality (almost exactly) can now be found inside Silo Cyc Obj.
//void HYMVectorObj::WriteData_SILO(DBfile *dbfile, int cycle,
//                                  char *mesh_name, float **mesh_coords) {
//    float *vec[ndims];
//    this->ReadVector_Binary(cycle,vec);
//    WriteVector_SILO(dbfile,this->varname,mesh_name,this->varnames,vec,
//                     mesh_coords[2],this->dims);
//    for(int m=0; m<ndims; m++)
//        delete [] vec[m];
//}
//
//============================================================================//
void HYMVectorObj::WriteData_ASCII(char *ascii_path, int cycle, double time, 
                                   float **mesh_coords) {
    char fname_ascii[1001];
    sprintf(fname_ascii,"%s_%03d.dat",this->ascii_name,cycle);
    float *vec[ndims];
    this->ReadVector_Binary(cycle,vec);
    WriteVector_ASCII(ascii_path,fname_ascii,this->vchar,this->dims,   
                      mesh_coords,time,vec,this->stopmsg);
    for(int m=0; m<ndims; m++)
        delete [] vec[m];
}

//============================================================================//
void HYMVectorObj::ReadVector_Binary(int cycle, float **vec) {
    this->PositionPointer_Binary(cycle);
    this->ReadVar_Binary(vec[0]);
    this->ReadVar_Binary(vec[1]);
    this->ReadVar_Binary(vec[2]);
}

//============================================================================//
//############################################################################//
//============================================================================//
