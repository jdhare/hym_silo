//============================================================================//
/*

Clayton Myers
SILO_CycleObj.cpp
Created:  06 February 2010
Modified: 06 February 2010

Class for SILO cycle objects..
                      
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
    this->data_objs = data_objs;
    this->stopmsg=stopmsg;
    
    this->time_str=NULL;
    this->stat_str="Clean";
    
    this->SetFlags(glob_data_flags);
    this->SetTime(); 
    if(halfcyl){
    	this->silodims = {dims[0],dims[1],dims[2]+1};
    }else{
    	this->silodims = {dims[0],dims[1],dims[2]+2*Nghost+1};
    }
    
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
    WriteMesh_SILO(dbfile,mesh_name,this->silodims,this->mesh_coords,
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

//Code previously in Silo write
//============================================================================//

void WriteMesh_SILO(DBfile *dbfile, char *mesh_name, int *dims, 
                    float **mesh_coords, int cycle, double time) {
    int i, j, k, n, Ntot;
    int Nq = dims[0], Nr = dims[1], Ns = dims[2];
    float *q = mesh_coords[0], *r = mesh_coords[1], *s = mesh_coords[2];
    float *s_ghost = NULL;
    

    //if halfcyl, just extend array. else ie. periodic, add ghost zones/
    if(half_cyl){
         Ns+=1;
         s_ghost=new float[Ns];
         //copy array across into new slightly alrger array
         for(k=0;k<Ns-1;k++){
                    s_ghost[k]=s[k];
         }
         //add final phi=pi point
         s_ghost[Ns-1]=s_ghost[Ns-2]+s_ghost[1];
    }else{
         AddGhostZones_Coord(s,s_ghost,Ns);               
    }
    Ntot = Nq*Nr*Ns;
    
    float *xg=NULL, *yg=NULL, *zg=NULL;
    xg = new float[Ntot];
    yg = new float[Ntot];
    zg = new float[Ntot];
    
    n = 0;
    for(k=0; k<Ns; k++) { 
        for(j=0; j<Nr; j++) {
            for(i=0; i<Nq; i++) {
                xg[n] = r[j]*cos(s_ghost[k]);
                yg[n] = r[j]*sin(s_ghost[k]);
                zg[n] = q[i];
                n++;
            } 
        }
    }
    
    delete [] s_ghost;
    
    // Create the arrays to write to the .silo database:
    int silodims[ndims] = {Nq,Nr,Ns};
    float *coords[ndims] = {(float*)xg, (float*)yg, (float*)zg};
    
    // Create an option list to save cycle and time values:
    DBoptlist *optlist = DBMakeOptlist(4);
    if(!half_cyl){
          int offset_low[ndims] = {0,0,Nghost};
          int offset_high[ndims] = {0,0,Nghost};
          DBAddOption(optlist, DBOPT_LO_OFFSET, offset_low);
          DBAddOption(optlist, DBOPT_HI_OFFSET, offset_high);                  
    }
    DBAddOption(optlist, DBOPT_DTIME, &time);
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);

    // DBAddOption(optlist, DBOPT_COORDSYS, DB_CYLINDRICAL);

    // Write the mesh to the .silo file:
    DBPutQuadmesh(dbfile, mesh_name, NULL, coords, silodims, ndims,
                  DB_FLOAT, DB_NONCOLLINEAR, optlist);
                  
    DBFreeOptlist(optlist);
    
    delete [] xg;
    delete [] yg;
    delete [] zg;
    
}

//============================================================================//
void WriteScalar_SILO(DBfile *dbfile, char *vname, char *mesh_name, float *var, 
                      int *dims) {
    // This function writes the data in var to the .silo database file
    int silodims[ndims] = {dims[0],dims[1],dims[2]+2*Nghost+1};
    float *silovar;
    if(halfcyl){
          AddFinalZone_Var(var,silovar,dims);              
    }else{
          // Add the SILO ghost zones to the data:
          AddGhostZones_Var(var,silovar,dims);
    }
    // Remove nonzero elements below a threshold (e.g. 1e-10):
    Set_Zeros(silovar,silodims[0]*silodims[1]*silodims[2]);
    // Write the data to the .silo file:
    DBPutQuadvar1(dbfile, vname, mesh_name, silovar, silodims, ndims, NULL, 0,
                  DB_FLOAT, DB_NODECENT, NULL);
    delete [] silovar;
}

//============================================================================//
void WriteVector_SILO(DBfile *dbfile, char *vname, char *mesh_name, 
                      char **varnames, float **vec, float *s, int* dims) {
    // This function writes the data in vec to the .silo database file
    int Ntot;
    int silodims[ndims] = {dims[0],dims[1],dims[2]+2*Nghost+1};
    float *silovec[ndims];
    if(halfcyl){
        AddFinalZone_Var(vec[0],silovec[0],dims);
        AddFinalZone_Var(vec[1],silovec[1],dims);  
        AddFinalZone_Var(vec[2],silovec[2],dims);  
    }else
        // Add the SILO ghost zones to the data:
        AddGhostZones_Var(vec[0],silovec[0],dims);
        AddGhostZones_Var(vec[1],silovec[1],dims);
        AddGhostZones_Var(vec[2],silovec[2],dims);
    }
    // Convert to Cartesian vector components:
    Cyl_to_Cart(silovec,s,dims);
    // Remove nonzero elements below a threshold (e.g. 1e-10):
    Ntot = silodims[0]*silodims[1]*silodims[2];
    Set_Zeros(silovec[0],Ntot);
    Set_Zeros(silovec[1],Ntot);
    Set_Zeros(silovec[2],Ntot);
    // Write the vector to the .silo file:
    DBPutQuadvar(dbfile, vname, mesh_name, ndims, varnames, silovec, silodims,
                 ndims, NULL, 0, DB_FLOAT, DB_NODECENT, NULL);    
    for(int m=0; m<ndims; m++) { delete [] silovec[m]; }
}

//============================================================================//
void AddGhostZones_Coord(float *s, float *&s_ghost, int &Ns) {
    int k, kshift;
    // Add the ghost zones in phi for mesh continuity:
    Ns += 2*Nghost + 1;
    s_ghost = new float[Ns];    
    kshift = Ns - (2*Nghost+1);
    for(k=Nghost; k<(Ns-Nghost-1); k++) {
        s_ghost[k] = s[k-Nghost];
        // Make periodicity continuous and add upper ghost zones:
        if(k <= (2*Nghost))
            s_ghost[k+kshift] = s_ghost[k];
        // Add lower ghost zones:
        if(k >= (Ns-(2*Nghost+1)))
            s_ghost[k-kshift] = s_ghost[k];
    }
}

//============================================================================//
void AddGhostZones_Var(float *var, float *&silovar, int *dims) {
    // Add the SILO ghost zones (in phi) to the variable
    int nshift, Ntot, n;
    int Nq = dims[0], Nr = dims[1], Ns = dims[2] + 2*Nghost + 1;
    
    Ntot = Nq*Nr*Ns;
    silovar = new float[Ntot];
    nshift = Ntot - (2*Nghost+1)*Nq*Nr;
    n = Nghost*Nq*Nr;
    for(int k=Nghost; k<(Ns-Nghost-1); k++) {
        for(int j=0; j<Nr; j++) {
            for(int i=0; i<Nq; i++) {
                silovar[n] = var[n-Nghost*Nq*Nr]; 
                // Continuous periodicity/ghost zones at top of mesh:
                if(k < (2*Nghost+1))
                    silovar[n+nshift] = silovar[n];
                // Ghost zones at the bottom of the mesh:
                if(k >= (Ns-(2*Nghost+1)))
                    silovar[n-nshift] = silovar[n];
                n++;                
            }
        }
    }
}

//============================================================================//
void AddFinalZone_Var(float *var, float *&silovar, int *dims) {
    // Add the SILO ghost zones (in phi) to the variable
    int nshift, Ntot, n;
    int Nq = dims[0], Nr = dims[1], Ns = dims[2] + 2*Nghost + 1;
    
    Ntot = Nq*Nr*Ns;
    silovar = new float[Ntot];
    nshift = Ntot - (2*Nghost+1)*Nq*Nr;
    n = Nghost*Nq*Nr;
    for(int k=Nghost; k<(Ns-Nghost-1); k++) {
        for(int j=0; j<Nr; j++) {
            for(int i=0; i<Nq; i++) {
                silovar[n] = var[n-Nghost*Nq*Nr]; 
                // Continuous periodicity/ghost zones at top of mesh:
                if(k < (2*Nghost+1))
                    silovar[n+nshift] = silovar[n];
                // Ghost zones at the bottom of the mesh:
                if(k >= (Ns-(2*Nghost+1)))
                    silovar[n-nshift] = silovar[n];
                n++;                
            }
        }
    }
}

//============================================================================//
void Cyl_to_Cart(float **vec, float *s, int *dims) {
    // Add functionality comment here
    int i, j, k, n, Ntot;
    int Nq = dims[0], Nr = dims[1], Ns = dims[2];
    float *vec_q = vec[0], *vec_r = vec[1], *vec_s = vec[2];
    float *vec_x=NULL, *vec_y=NULL, *s_ghost=NULL;
    
    AddGhostZones_Coord(s,s_ghost,Ns);
    Ntot = Nq*Nr*Ns;
    vec_x = new float[Ntot];
    vec_y = new float[Ntot];
    
    // Convert (r,phi) components to (x,y) components
    n = 0;
    for(k=0; k<Ns; k++) {
        for(j=0; j<Nr; j++) {
            for(i=0; i<Nq; i++) {
                vec_x[n] = vec_r[n]*cos(s_ghost[k]) - vec_s[n]*sin(s_ghost[k]);
                vec_y[n] = vec_r[n]*sin(s_ghost[k]) + vec_s[n]*cos(s_ghost[k]);
                n++;
            }
        }
    }
    //set exclude origin in HYM_SILO.hpp
    if(exclude_origin){
        // Interpolate on axis at r=0
        for(i=0; i<Nq; i++) {
            int nsum=0;
            float vec_x_sum=0., vec_y_sum=0.;
            for(k=Nghost; k<(Ns-Nghost-1); k++) {
                n = fn(i,1,k,Nq,Nr);
                vec_x_sum += vec_x[n];
                vec_y_sum += vec_y[n];
                nsum++;
            }
            for(k=0; k<Ns; k++) {
                n = fn(i,0,k,Nq,Nr); 
                vec_x[n] = vec_x_sum/nsum;
                vec_y[n] = vec_y_sum/nsum;
            }
        }
    }
    
    // Replace and sort the vector components into the (x,y,z) arrangement
    delete [] vec_r; delete [] vec_s; delete [] s_ghost;
    vec[0] = vec_x;
    vec[1] = vec_y;
    vec[2] = vec_q;
}

//============================================================================//
void Set_Zeros(float *var, int Ntot) {
    float dum;
    for(int n=0; n<Ntot; n++) {
        if(var[n] < 0)
            dum = -var[n];
        else
            dum = var[n];
        if(dum < zero_threshold)
            var[n] = 0.0;
    }
}

//============================================================================//
//============================================================================//
