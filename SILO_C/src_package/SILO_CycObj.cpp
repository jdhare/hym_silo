//============================================================================//
/*

Clayton Myers
SILO_CycleObj.cpp
Created:  06 February 2010
Modified: 06 February 2010

Class for SILO cycle objects.
                      
*/
//============================================================================//
//============================================================================//

#include <HYM_SILO.hpp>
#include <HYM_DataObj.hpp>
#include <SILO_Write.hpp>
#include <SILO_CycObj.hpp>

//============================================================================//
//============================================================================//

char *silo_name = "HYM";
char *mesh_name = "HYM_mesh";

//============================================================================//
//============================================================================//
SILO_CycObj::SILO_CycObj(int cycle, float **mesh_coords, int *dims, 
                         HYMDataObj **data_objs, bool *glob_data_flags, 
                         char *stopmsg) {
    this->cycle = cycle;
    this->mesh_coords = mesh_coords;
    this->dims = dims;
    this->data_objs = data_objs;
    this->stopmsg=stopmsg;
    
    this->time_str=NULL;
    this->stat_str="Clean";
    
    this->SetFlags(glob_data_flags);
    this->SetTime();    
}

//============================================================================//
SILO_CycObj::~SILO_CycObj(void) { 
    return;
}


//============================================================================//
void SILO_CycObj::SetFlags(bool *glob_data_flags) {
    this->write_flag  = false;
    this->report_flag = false;

    for(int m=0; m<nvars; m++) {
        if(this->data_objs[m] != NULL) { 
            int Ncyc_m = this->data_objs[m]->Ncyc;
	    if(this->cycle <= Ncyc_m)
                this->mask_flags[m] = this->data_objs[m]->cycle_mask[this->cycle-1];
            else 
                this->mask_flags[m] = false;
        }
        else
            this->mask_flags[m] = false;
        if(this->mask_flags[m])
            this->write_flag = true;
        if(this->mask_flags[m] != glob_data_flags[m]) {
            this->report_flag = true;
            this->stat_str = "Data";
        }
    }
    if(!this->write_flag)
        this->time_str = "NULL";
}

//============================================================================//
void SILO_CycObj::SetTime(void) {
    this-> time = 0.0;
    bool found_first_time_value = false;
    for(int m=0; m<nvars; m++) {
        if(this->mask_flags[m]) {
            if(!found_first_time_value) {
                time = this->data_objs[m]->times[cycle-1];
                found_first_time_value = true;
            }
            else if(time != this->data_objs[m]->times[cycle-1]) {
                time = -1.0;
                this->report_flag = true;
                this->stat_str = "Time";
                this->time_str = "NULL";
            }
        }
    }
}

//============================================================================//
void SILO_CycObj::Write_SILO(char *silo_path) {
    // Write out the SILO database filename:
    char full_name[1001];
    sprintf(full_name,"%s%s_%0.3d.silo",silo_path,silo_name,this->cycle);

    // Bail out if this SILO database should not be written:
    if(!write_flag) {
        cout << "      Ignored: " << full_name << "\n";
        return;
    }
    
    // Open the .silo database file:
    DBfile *dbfile = DBCreate(full_name,DB_CLOBBER,DB_LOCAL,"data",DB_PDB);   
    // Write the mesh to the .silo database:
    WriteMesh_SILO(dbfile,mesh_name,this->dims,this->mesh_coords,
                   this->cycle,this->time);    
    // Write the data to the .silo database:
    for(int m=0; m<nvars; m++) {
        if(this->mask_flags[m]) {
            this->data_objs[m]->WriteData_SILO(dbfile,this->cycle,mesh_name,
                                               this->mesh_coords);
        }
    }
    // Close the completed .silo database:
    DBClose(dbfile);
    cout << "      Output:  " << full_name << "\n";    
}

//============================================================================//
void SILO_CycObj::Write_ASCII(char *ascii_path) {
    VerifyPath(ascii_path,stopmsg);
    for(int m=0; m<nvars; m++) {
        if(this-mask_flags[m]) {
            this->data_objs[m]->WriteData_ASCII(ascii_path,this->cycle,
                                                this->time,this->mesh_coords);
        }
    }                 
}

//============================================================================//
//============================================================================//
