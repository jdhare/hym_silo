//============================================================================//
/*

Clayton Myers
SILO_Read.hpp
Created:  17 September 2009
Modified: 17 September 2009

Header file for SILO read functions.

*/
//============================================================================//
//============================================================================//

void OpenFile_SILO(char*,char*,DBfile*&,char*);
void ReadScalar_SILO(char*,char*,char*,float*&,int*,char*);
void ReadVector_SILO(char*,char*,char*,float**,int*,char*);
double ReadTime_SILO(char*,char*,char*,char*);
int Get_Ncyc(char*,char*);
void Get_Mesh_Dims(char*,char*,char*,int*,char*);

//============================================================================//
//============================================================================//
