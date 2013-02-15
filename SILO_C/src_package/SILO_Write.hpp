//============================================================================//
/*

Clayton Myers
SILO_Write.hpp
Created:  17 September 2009
Modified: 17 September 2009

Header file for SILO write functions.

*/
//============================================================================//
//============================================================================//

void WriteMesh_SILO(DBfile*,char*,int*,float**,int,double);
void WriteScalar_SILO(DBfile*,char*,char*,float*,int*);
void WriteVector_SILO(DBfile*,char*,char*,char**,float**,float*,int*);

//============================================================================//
//============================================================================//
