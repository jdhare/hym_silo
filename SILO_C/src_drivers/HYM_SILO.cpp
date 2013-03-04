//============================================================================//
/*

Clayton Myers
HYM_SILO.cpp
Created:  18 June 2009
Modified: 06 February 2010

This is the parent routine to the HYM-to-SILO conversion routine.  It accesses 
the raw binary output from HYM and converts one or all of the output cycles to
SILO databases.

Command line arguments:
    (1) data_path  -- Path to the location of the raw binary data files.
                      Must include leading and trailing slashes.
    (2) silo_path  -- Path to the destination for the .silo database(s).
                      Must include leading and trailing slashes.
    (3) cycle      -- Cycle number to be converted in the range
                      [1,Ncyc].  If cycle = 0, the entire range is converted.
    (4) data_flags -- String of 1's and 0's to mask which variables are added
                      to the converted .silo database(s).
                      
Present data_flags key: (Length = 5)
    #1 = p (pressure)
    #2 = n (density)
    #3 = B (magnetic field)
    #4 = v (fluid velocity)
    #5 = J (current density)
                      
*/
//============================================================================//
//============================================================================//

#include <HYM_SILO.hpp>
#include <HYM_DataObj.hpp>
#include <SILO_CycObj.hpp>

//============================================================================//
//============================================================================//
void ReadArgs(int,char**,char*&,char*&,int&,bool*);
void ReadStatData(char*,int&,int*);
void Unify_Ncyc(HYMDataObj**,int&);
void Write_Report(int Ncyc,SILO_CycObj**,char*);
void CleanUp(int,float**,HYMDataObj**,SILO_CycObj**);

char *fname_mesh = "hgrid.d";
char *fname_stat = "hstat.d";
char *stopmsg = "Stopping SILO File Construction.";

//============================================================================//
int main(int argc, char *argv[]) {
    int cycle, Ncyc, dims[ndims];
    float *mesh_coords[ndims] = {NULL,NULL,NULL};
    char *data_path=NULL, *silo_path=NULL;
    bool data_flags[nvars], report_flag;
    
    //------------------------------------------------------------------------//
    // Process the basic run parameters and read in the mesh:
    ReadArgs(argc,argv,data_path,silo_path,cycle,data_flags);
    ReadStatData(data_path,Ncyc,dims);
    HYMDataObj::ReadMesh_Binary(data_path,fname_mesh,dims,mesh_coords,stopmsg);
    cout << "Code made it here."<<endl;
    
    // Define and initialize the HYM data objects:
    HYMDataObj *data_objs[nvars];
    for(int m=0; m<nvars; m++)
        data_objs[m] = NULL;
    
    // Use data_flags to determine which data values should be included:
    if(data_flags[0])
        data_objs[0] = new HYMScalarObj('p',data_path,"h3ds.d",dims,stopmsg);
    if(data_flags[1])
        data_objs[1] = new HYMScalarObj('n',data_path,"h3ds_ff.d",dims,stopmsg);
    if(data_flags[2])
        data_objs[2] = new HYMVectorObj('B',data_path,"h3db.d",dims,stopmsg);
    if(data_flags[3])
        data_objs[3] = new HYMVectorObj('v',data_path,"h3dv.d",dims,stopmsg);
    if(data_flags[4])
        data_objs[4] = new HYMVectorObj('J',data_path,"h3dj.d",dims,stopmsg);

    // Unify the Ncyc value and the time series:
    Unify_Ncyc(data_objs,Ncyc);
    
    // Define and initialize the SILO cycle objects:
    SILO_CycObj *cyc_objs[Ncyc];
    for(int m=0; m<Ncyc; m++)
        cyc_objs[m] = NULL;
    
    //------------------------------------------------------------------------//
    // Construct the SILO cycle objects:
    if(cycle == 0) {
        for(int m=0; m<Ncyc; m++)
            cyc_objs[m] = new SILO_CycObj(m+1,mesh_coords,dims,data_objs,
                                          data_flags,stopmsg);
    }
    else if(cycle >= 1 && cycle <= Ncyc)
        cyc_objs[cycle-1] = new SILO_CycObj(cycle,mesh_coords,dims,data_objs,
                                            data_flags,stopmsg);
    else {
        CleanUp(Ncyc,mesh_coords,data_objs,cyc_objs);
        char message[1001];
        sprintf(message,"      %s%d%s\n      %s%d%s\n      %s",
                "The requested cycle number is not valid (cycle = ",cycle,").",
                "The cycle number must be zero (0) or in the range [1,",Ncyc,
                "].",stopmsg);
        StopExecution(message);
    }

    // Write the SILO databases:
    report_flag = false;  
    for(int m=0; m<Ncyc; m++) {
        if(cyc_objs[m] != NULL) {
            cyc_objs[m]->Write_SILO(silo_path);
            if(cyc_objs[m]->report_flag)
                report_flag = true;   
        }
    }
    
    // Write a report if abnormalities exist:
    if(report_flag)
        Write_Report(Ncyc,cyc_objs,silo_path);    

    // Clean up the allocated arrays:
    CleanUp(Ncyc,mesh_coords,data_objs,cyc_objs);
}

//============================================================================//
//============================================================================//
void ReadArgs(int argc, char **argv, char *&data_path, char *&silo_path, 
              int &cycle, bool *data_flags) {
    char message[1001];
              
    // Count the initial command line arguments
    if(!(argc == 5)) {
        sprintf(message,"      %s%s\n      %s","An improper number of ",
                "command line arguments was found.",stopmsg);
        StopExecution(message);
    }
        
    // Distribute the command line arguments:
    data_path = argv[1];
    silo_path = argv[2];
    VerifyPath(data_path,stopmsg);
    VerifyPath(silo_path,stopmsg);
    ConvertToInt(argv[3],cycle,stopmsg);
    
    // Determine if the final command line argument has the proper length:
    if(argv[4][nvars] != '\0') {
        sprintf(message,"      %s%s\n         %s %d; %s %s\n      %s",
                "The data_flags command line argument does not have the ", 
                "proper length:", "Required Length:",nvars,
                "Received Argument:",argv[4],stopmsg);
        StopExecution(message);
    }
    // Convert the final argument to a boolean array of length nvars:
    for(int m=0; m<nvars; m++) {
        if(argv[4][m] == '1')
            data_flags[m] = true;
        else if(argv[4][m] == '0')
            data_flags[m] = false;
        else {
            sprintf(message,"      %s\n  %s%d;  %s%c.\n        %s\n      %s",
                "An improper data flag was found:",
                "      Flag Slot = ",m,"Flag = ", argv[4][m],
                "Data flags must be either 0 or 1.",stopmsg);
            StopExecution(message);
        }
    }
}

//============================================================================//
void ReadStatData(char *data_path, int &Ncyc, int *dims) {
    // Gets the number of cycles (Ncyc) and mesh dimensions (dims) of the run.
    char str[1001];
    ifstream file;
    
    // Verify the status file:
    if(!VerifyInputFile(data_path,fname_stat,"binary",stopmsg)) {
        char message[1001];
        sprintf(message,"      %s\n      %s%s\n      %s%s\n      %s",
                "The status file was not found:","    Data Path:    ",
                data_path,"    Missing File: ",fname_stat, stopmsg);
        StopExecution(message);
    }
    
    // Open the status file:
    OpenInputFile(file,data_path,fname_stat,"ascii",stopmsg);
    
    // Scroll:
    for(int m=0; m<4; m++) { file >> str; }
    // Save the dimensions:
    for(int m=0; m<ndims; m++) { file >> str >> dims[m]; }
    // Scroll:
    for(int m=0; m<17; m++) { file >> str; }
    // Save the number of cycles:
    file >> Ncyc;
    
    file.close();
    
    sprintf(str,"%-10.10s",fname_stat);
    cout << "      Reading: " << str;
    cout << " Ncyc = " << Ncyc << endl;
    
    // Add back the ghost zones in r near the origin:
    dims[1] += 2;
}

//============================================================================//
void Unify_Ncyc(HYMDataObj **data_objs, int &Ncyc) {
    // This function compares the Ncyc values from each of the HYM data objects
    // and sends warnings if they are not the same.  If more data exists than
    // is indicated by hstat.d, the global Ncyc value is increased accordingly.
    int Ncyc_max = 0;
    for(int m=0; m<nvars; m++) {
        if(data_objs[m] != NULL && data_objs[m]->Ncyc != Ncyc) {
            cout << "\n      Warning: Ncyc for the " << data_objs[m]->vchar;
            cout << " data object\n";
            cout << "               differs from the hstat.d Ncyc value.\n\n"; 
        }
        if(data_objs[m] != NULL && data_objs[m]->Ncyc > Ncyc_max)
            Ncyc_max = data_objs[m]->Ncyc;
    }
    if(Ncyc_max == 0) {
        cout << "\n      Warning: No data was found to convert to the SILO ";
        cout << "format.\n";
        cout << "                 Only a mesh will be written to each SILO ";
        cout << "database.\n\n";
    }
    else if(Ncyc_max != Ncyc) {
        Ncyc = Ncyc_max;
        cout << "\n      Warning: The global Ncyc value has been changed to ";
        cout << "Ncyc = " << Ncyc << ".\n\n";
    }
}

//============================================================================//
void Write_Report(int Ncyc, SILO_CycObj **cyc_objs, char *silo_path) {
    cout << "\n      Abnormal conditions were detected during SILO conversion.";
    cout << "\n      A report of the output of this run has been generated.";
    cout << "\n      It has been saved as \"SILO_Report.n\" ";
    cout << "in the SILO output directory.\n\n";
    
    ofstream file;
    char outstr[1001];    
    OpenOutputFile(file,silo_path,"SILO_Report.n",stopmsg);
    file << "Cycle      Time    Out?     pnBvJ     Fault\n";
    for(int m=0; m<Ncyc; m++) {
        if(cyc_objs[m] != NULL) {
            sprintf(outstr,"  %03d",cyc_objs[m]->cycle);
            file << outstr;
            if(cyc_objs[m]->time_str == NULL)
                sprintf(outstr,"%10.1f",cyc_objs[m]->time);
            else
                sprintf(outstr,"%10s",cyc_objs[m]->time_str);
            file << outstr;
            if(cyc_objs[m]->write_flag)
                sprintf(outstr,"%8s","Yes");
            else
                sprintf(outstr,"%8s","No");
            file << outstr << "     ";
            for(int mm=0; mm<nvars; mm++)
                file << cyc_objs[m]->mask_flags[mm];
            sprintf(outstr,"%10s",cyc_objs[m]->stat_str);
            file << outstr << endl;
        }
    }
    file.close();
}

//============================================================================//
void CleanUp(int Ncyc, float **mcoords, HYMDataObj **dobj, SILO_CycObj **cobj) {
    for(int m=0; m<ndims; m++)
        delete [] mcoords[m];
    for(int m=0; m<nvars; m++) {
        if(dobj[m] != NULL)
            delete dobj[m];
    }
    for(int m=0; m<Ncyc; m++) {
        if(cobj[m] != NULL)
            delete cobj[m];
    }
}

//============================================================================//
//============================================================================//
