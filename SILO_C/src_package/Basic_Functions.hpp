//============================================================================//
/*

Clayton Myers
Basic_Functions.hpp
Created:  17 September 2009
Modified: 11 December 2009

Header file for basic functions.

*/
//============================================================================//
//============================================================================//

void OpenInputFile(ifstream&,char*,char*,char*,char*);
bool VerifyInputFile(char*,char*,char*,char*);
void OpenOutputFile(ofstream&,char*,char*,char*);
void OpenFile_RW(fstream&,char*,char*,char*);
void ClearFile(char*,char*,char*);
void VerifyPath(char*,char*);
void ConvertToInt(char*,int&,char*);
void ConvertToFloat(char*,float&,char*);
void ConvertToBool(char*,bool&,char*);
void StopExecution(char*);
void PrintTime(time_t);
void Print_Dims(int*);
void Print_Dims(int*,char*);
void Construct_Phi(float*&,int*);
int fn(int,int,int,int,int);

//============================================================================//
void listdir(char*,char**&,int&,char*);

//============================================================================//
//============================================================================//
