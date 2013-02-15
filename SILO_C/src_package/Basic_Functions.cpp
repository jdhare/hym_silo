//============================================================================//
/*

Clayton Myers
Basic_Functions.cpp
Created:  12 August 2009
Modified: 06 February 2010

Support functions for performing basic operations and error checking.

*/
//============================================================================//
//============================================================================//

#include <HYM_SILO.hpp>

void OpenFile_Safe(ifstream&,char*,char*,char*);

//============================================================================//
void OpenInputFile(ifstream &file, char *fpath, char *fname, char *format,
                   char *stopmsg) {                   
    char path[1001], message[1001];
    strcpy(path,fpath);
    strcat(path,fname);
    OpenFile_Safe(file,path,format,stopmsg);
     
    if(file.fail()) {
        sprintf(message,"      %s\n      %s%s%s\n      %s",
        "Error in function OpenInputFile.", 
        "The file \"",path,"\" could not be found.",stopmsg);
        StopExecution(message);        
    }
}

//============================================================================//
bool VerifyInputFile(char *fpath, char *fname, char *format, char *stopmsg) {
    ifstream file;
    char path[1001];
    strcpy(path,fpath);
    strcat(path,fname);
    OpenFile_Safe(file,path,format,stopmsg);
        
    bool verified = !file.fail();
    file.close();
    return verified;
}

//============================================================================//
void OpenFile_Safe(ifstream &file, char *path, char *format, char *stopmsg) {
    char message[1001];
    if(strcmp(format,"ascii") == 0)
        file.open(path);
    else if(strcmp(format,"binary") == 0)
        file.open(path,ios::in|ios::binary);
    else {
        sprintf(message,"      %s\n  %s%s%s\n  %s%s%s\n      %s",
                "Error in function OpenInputFile.", 
                "The file format string \"",format,"\" is not valid.",
                "The file \"",path,"\" could not be opened.",stopmsg);
        StopExecution(message);   
    }
}

//============================================================================//
void OpenOutputFile(ofstream &file, char *fpath, char *fname, char *stopmsg) {
    char path[1001], message[1001];
    strcpy(path,fpath);
    strcat(path,fname);
    file.open(path,ios::out);

    if(file.fail()) {
        sprintf(message,"      %s\n  %s%s%s\n      %s",
        "Error in function OpenOutputFile.", 
        "The file \"",path,"\" could not be found.",stopmsg);
        StopExecution(message);        
    }
}

//============================================================================//
void OpenFile_RW(fstream &file, char *fpath, char *fname, char *stopmsg) {
    char path[1001], message[1001];
    strcpy(path,fpath);
    strcat(path,fname);
    file.open(path);

    if(file.fail()) {
        sprintf(message,"      %s\n  %s%s%s\n      %s",
        "Error in function OpenFile_RW.", 
        "The file \"",path,"\" could not be found.",stopmsg);
        StopExecution(message);        
    }
}

//============================================================================//
void ClearFile(char *fpath, char *fname, char *stopmsg) {
    char path[1001], message[1001];
    ofstream file;
    strcpy(path,fpath);
    strcat(path,fname);
    file.open(path,ios::out);
    
    if(file.fail()) {
        sprintf(message,"      %s\n  %s%s%s\n      %s",
        "Error in function ClearFile.", 
        "The file \"",path,"\" could not be found.",stopmsg);
        StopExecution(message);        
    }
    
    file.close();
}

//============================================================================//
void VerifyPath(char *path, char *stop_message) {
    // Verifies the existence of the directory given by the path argument
    char orig_dir[1001];
    getcwd(orig_dir,1001);
    if(chdir(path) != 0) {
        char message[1001];
        sprintf(message,"      The directory \"%s\" %s\n      %s",path,
                "could not be found.",stop_message);
        StopExecution(message);
    }
    chdir(orig_dir);
}

//============================================================================//
void ConvertToInt(char *arg, int &argi, char *stopmsg) {
    // Converts the string argument arg to an integer with error checking
    // It is assumed that there are no leading or trailing white spaces in arg
    stringstream ssout(arg);
    if(!(ssout >> argi)) {
        char message[1001];
        sprintf(message,"      The argument \"%s\" %s\n      %s",arg,
                "could not be converted to an integer.",stopmsg);
        StopExecution(message);
    }
}

//============================================================================//
void ConvertToFloat(char *arg, float &argf, char *stopmsg) {
    // Converts the string argument arg to a float with error checking
    stringstream ssout(arg);
    if (!(ssout >> argf)) {
        char message[1001];
        sprintf(message,"      The argument \"%s\" %s\n      %s",arg,
                "could not be converted to a float.",stopmsg);
        StopExecution(message);
    }
}

//============================================================================//
void ConvertToBool(char *arg, bool &argf, char *stopmsg) {
    // Converts the string argument arg to a float with error checking
    stringstream ssout(arg);
    if (!(ssout >> argf)) {
        char message[1001];
        sprintf(message,"      The argument \"%s\" %s\n      %s",arg,
                "could not be converted to a bool.",stopmsg);
        StopExecution(message);
    }
}

//============================================================================//
void StopExecution(char *message) {
    printf("\n%s\n\n",message);
    // printf("%s\n",message);
    exit(0);
}

//============================================================================//
void PrintTime(time_t start_time) {
    // Prints the execution time in whole seconds.
    int time, min, sec;
    
    time_t end_time = std::time(NULL);
    
    time = (int)(end_time - start_time);
    min = (int)floor(time/60.);
    sec = (int)(time-60*min);
    
    printf("\n  Execution Time: %02d:%02d\n\n",min,sec);
}

//============================================================================//
void Print_Dims(int *dims) {
    printf("\n");
    for(int m=0; m<ndims; m++)
        printf("  dims[%d] = %3d\n",m,dims[m]);
    printf("\n");
}

//============================================================================//
void Print_Dims(int *dims, char *buffer) {
    printf("\n");
    for(int m=0; m<ndims; m++) {
        printf(buffer);
        printf("  dims[%d] = %3d\n",m,dims[m]);
    }
    printf("\n");
}

//============================================================================//
void Construct_Phi(float *&phi, int *dims) {
    phi = new float[dims[2]];
    float delta = 2.0*pi/dims[2];
    for(int k=0; k<dims[2]; k++)
        phi[k] = k*delta;        
}

//============================================================================//
int fn(int i, int j, int k ,int Nz, int Nr) {
    // Transforms the 3D indices (i,j,k) to the linear index n
    return k*Nr*Nz + j*Nz + i; 
}

//============================================================================//
//============================================================================//

#include <dirent.h>
#include <cstdio>

void listdir(char *path, char **&dirlist, int &nlist, char *stopmsg) {
    DIR *pdir = NULL;
    struct dirent *pent = NULL;
    
    VerifyPath(path,stopmsg);
    
    // Load the directory and count the entries:
    nlist = 0;
    pdir = opendir(path);
    if(pdir == NULL)
        StopExecution("ERROR! pdir could not be initialized correctly");
    while(pent = readdir(pdir)) {
        if(pent == NULL)
            StopExecution("ERROR! pent could not be initialized correctly");
        nlist++;
    }
    closedir(pdir);
    
    // Allocate the memory for the entry strings:
    dirlist = new char*[nlist];
    
    // Load the directory again and save the entries:
    nlist = 0;
    pent = NULL;
    pdir = opendir(path);
    if(pdir == NULL)
        StopExecution("ERROR! pdir could not be initialized correctly");
    while(pent = readdir(pdir)) {
        if(pent == NULL)
            StopExecution("ERROR! pent could not be initialized correctly");
        dirlist[nlist] = pent->d_name;
        nlist++;
    }
    closedir(pdir);
}

//============================================================================//
//============================================================================//

