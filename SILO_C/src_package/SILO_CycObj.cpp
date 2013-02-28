//============================================================================//
/*

Clayton Myers
SILO_CycleObj.cpp
Created:  06 February 2010
Modified: 06 February 2010
 * Modified: 28 February 2012. Rewrite to include Silo_write functions

Class for SILO cycle objects.
                      
*/
//============================================================================//
//============================================================================//

#include <HYM_SILO.hpp>
#include <HYM_DataObj.hpp>
#include <SILO_CycObj.hpp>

//============================================================================//
//============================================================================//

char *silo_name = "HYM";
char *mesh_name = "HYM_mesh";
const float zero_threshold = 1.0E-7; // For the Set_Zero function

//============================================================================//
//============================================================================//
SILO_CycObj::SILO_CycObj(int cycle, float **mesh_coords, int *dims, 
                         HYMDataObj **data_objs, bool *glob_data_flags, 
                         char *stopmsg) {
    this->cycle = cycle;
    this->data_objs = data_objs;
    this->stopmsg=stopmsg;
    this->silo_mesh=new float*[3];
    
    this->time_str=NULL;
    this->stat_str="Clean";
    
    this->SetFlags(glob_data_flags);
    this->SetTime();
    //half_cyl defined in HYM_SILO.hpp - at the moment, we either half a half cylinder or a full, periodic cylinder
    if(half_cyl){
        this->AddFinalZone(dims, mesh_coords);
    }else{
        this->AddGhostZones(dims, mesh_coords);
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
void SILO_CycObj::AddFinalZone(int *dims, float **mesh_coords){
    //increase the dimensions by one
    this->silodims=dims;
    this->silodims[2]+=1;
    //copy the unmodified coordinates over
    this->silo_mesh[0]=mesh_coords[0];
    this->silo_mesh[1]=mesh_coords[1];
    //add the final mesh coordinate in phi
    this->silo_mesh[2]=new float[silodims[2]];
    for(int k=0; k<silodims[2]; k++){
        //very explicit - make the last coord one increment greater than the previous coord
        if(k==(silodims[2]-1)){
            silo_mesh[2][k]=mesh_coords[2][k-1]+(mesh_coords[2][1]-mesh_coords[2][0]);
        }else{
            silo_mesh[2][k]=mesh_coords[2][k];
        }
    }
}
//============================================================================//
void SILO_CycObj::AddGhostZones(int *dims, float **mesh_coords){
    //increase the dimensions by one
    this->silodims=dims;
    this->silodims[2]+=2*Nghost+1;
    //copy the unmodified coordinates over
    this->silo_mesh[0]=mesh_coords[0];
    this->silo_mesh[1]=mesh_coords[1];
    //add the final mesh coordinate in phi
    this->silo_mesh[2]=new float[silodims[2]];
    int kshift = silodims[2] - (2*Nghost+1);
    for(int k=Nghost;k<silodims[2]-Nghost-1;k++){
        silo_mesh[2][k]=mesh_coords[2][k];
        if(k<=2*Nghost)
            silo_mesh[2][k+kshift]=mesh_coords[2][k];
        if(k>=kshift)
            silo_mesh[2][k-kshift]=mesh_coords[2][k];
        
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
    this->WriteMesh_SILO(dbfile);    
    // Write the data to the .silo database:
    for(int m=0; m<nvars; m++) {
        if(this->mask_flags[m]) {
            if(data_objs[m]->id=='S'){
                //scalar stuff
                float* vars=NULL;
                data_objs[m]->GetData_SILO(this->cycle,&vars);
                this->WriteScalar_SILO(dbfile,data_objs[m]->varname,*vars);
                delete vars;
            }
            if(data_objs[m]->id=='V'){
                //vector stuff
                float** vecs=NULL;
                data_objs[m]->GetData_SILO(this->cycle,&vecs);
                this->WriteVector_SILO(dbfile,data_objs[m]->varname,data_objs[m]->varnames,vecs);
                for(int m=0; m<ndims; m++)
                    delete [] vecs[m];
            }
        }
    }
    // Close the completed .silo database:
    DBClose(dbfile);
    cout << "      Output:  " << full_name << "\n";    
}
//============================================================================//
void SILO_CycObj::Write_ASCII(char *ascii_path, float **mesh_coords) {
    VerifyPath(ascii_path,stopmsg);
    for(int m=0; m<nvars; m++) {
        if(this-mask_flags[m]) {
            this->data_objs[m]->WriteData_ASCII(ascii_path,this->cycle,
                                                this->time,**mesh_coords);
        }
    }                 
}

//============================================================================//
void SILO_CycObj::WriteMesh_SILO(DBfile *dbfile) {
    int i, j, k, n, Ntot;
    int Nq = silodims[0], Nr = silodims[1], Ns = silodims[2];
    float *q = silo_mesh[0], *r = silo_mesh[1], *s = silo_mesh[2];

    Ntot = Nq*Nr*Ns;
    
    float *xg=NULL, *yg=NULL, *zg=NULL;
    xg = new float[Ntot];
    yg = new float[Ntot];
    zg = new float[Ntot];
    
    //convert to Cartesian. n is ugly, we often use fn() yet, but I haven't rewritten this
    n = 0;
    for(k=0; k<Ns; k++) { 
        for(j=0; j<Nr; j++) {
            for(i=0; i<Nq; i++) {
                xg[n] = r[j]*cos(s[k]);
                yg[n] = r[j]*sin(s[k]);
                zg[n] = q[i];
                n++;
            } 
        }
    }
        
    // Create the arrays to write to the .silo database: (ndims is a global)
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
void SILO_CycObj::WriteScalar_SILO(DBfile *dbfile, char *vname, float *var){
    // This function writes the data in var to the .silo database file
    float *silovar;
    if(half_cyl){
          this->AddFinalZone_Var(var,silovar);              
    }else{
          // Add the SILO ghost zones to the data:
          this->AddGhostZones_Var(var,silovar);
    }
    // Remove nonzero elements below a threshold (e.g. 1e-10):
    Set_Zeros(silovar,silodims[0]*silodims[1]*silodims[2]);
    // Write the data to the .silo file:
    DBPutQuadvar1(dbfile, vname, mesh_name, silovar, silodims, ndims, NULL, 0,
                  DB_FLOAT, DB_NODECENT, NULL);
    delete [] silovar;
}

//============================================================================//
void SILO_CycObj::WriteVector_SILO(DBfile *dbfile, char *vname, 
                      char **varnames, float **vec) {
    // This function writes the data in vec to the .silo database file
    int Ntot;
    float *silovec[ndims];
    if(half_cyl){
        this->AddFinalZone_Var(vec[0],silovec[0]);
        this->AddFinalZone_Var(vec[1],silovec[1]);  
        this->AddFinalZone_Var(vec[2],silovec[2]);  
    }else{
        // Add the SILO ghost zones to the data:
        this->AddGhostZones_Var(vec[0],silovec[0]);
        this->AddGhostZones_Var(vec[1],silovec[1]);
        this->AddGhostZones_Var(vec[2],silovec[2]);
    }
    // Convert to Cartesian vector components:
    this->Cyl_to_Cart(silovec);
    // Remove nonzero elements below a threshold (e.g. 1e-10):
    Ntot = silodims[0]*silodims[1]*silodims[2];
    this->Set_Zeros(silovec[0],Ntot);
    this->Set_Zeros(silovec[1],Ntot);
    this->Set_Zeros(silovec[2],Ntot);
    // Write the vector to the .silo file:
    DBPutQuadvar(dbfile, vname, mesh_name, ndims, varnames, silovec, this->silodims,
                 ndims, NULL, 0, DB_FLOAT, DB_NODECENT, NULL);    
    for(int m=0; m<ndims; m++) { delete [] silovec[m]; }
}

//============================================================================//
void SILO_CycObj::AddGhostZones_Var(float *var, float *&silovar) {
    // Add the SILO ghost zones (in phi) to the variable
    int nshift, Ntot, n;
    int Nq = silodims[0], Nr = silodims[1], Ns = silodims[2];
   
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
void SILO_CycObj::AddFinalZone_Var(float *var, float *&silovar) {
    // Add the SILO ghost zones (in phi) to the variable
    int nshift, Ntot, n;
    int Nq = silodims[0], Nr = silodims[1], Ns = silodims[2];
    
    Ntot = Nq*Nr*Ns;
    silovar = new float[Ntot];
    nshift = Ntot - Nq*Nr;
    n = 0;
    for(int k=0; k<Ns-1; k++) {
        for(int j=0; j<Nr; j++) {
            for(int i=0; i<Nq; i++) {
                silovar[n] = var[n]; 
                // Copy variables from phi=0 to phi=pi
                if(k ==0)
                    silovar[n+nshift] = silovar[n];
                n++;                
            }
        }
    }
}

//============================================================================//
void SILO_CycObj::Cyl_to_Cart(float **vec) {
    // Add functionality comment here
    int i, j, k, n, Ntot;
    int Nq = silodims[0], Nr = silodims[1], Ns = silodims[2];
    float *vec_q = vec[0], *vec_r = vec[1], *vec_s = vec[2];
    float *vec_x=NULL, *vec_y=NULL;
    
    Ntot = Nq*Nr*Ns;
    vec_x = new float[Ntot];
    vec_y = new float[Ntot];
    
    // Convert (r,phi) components to (x,y) components
    n = 0;
    for(k=0; k<Ns; k++) {
        for(j=0; j<Nr; j++) {
            for(i=0; i<Nq; i++) {
                vec_x[n] = vec_r[n]*cos(silo_mesh[2][k]) - vec_s[n]*sin(silo_mesh[2][k]);
                vec_y[n] = vec_r[n]*sin(silo_mesh[2][k]) + vec_s[n]*cos(silo_mesh[2][k]);
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
    delete [] vec_r; delete [] vec_s;
    vec[0] = vec_x;
    vec[1] = vec_y;
    vec[2] = vec_q;
}

//============================================================================//
void SILO_CycObj::Set_Zeros(float *var, int Ntot) {
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
//============================================================================//
