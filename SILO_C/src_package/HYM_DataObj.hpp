//============================================================================//
/*

Clayton Myers
HYM_DataObj.hpp
Created:  15 September 2009
Modified: 06 February 2010

Header file for HYM scalar and vector data objects.

*/
//============================================================================//
//============================================================================//
class HYMDataObj {
    public:
        char vchar;          // Single character name of the variable
        int Ncyc;            // Number of cycles in the source data
        bool *cycle_mask;    // Mask of cycles where the data exists
        double *times;       // Vector with the simulation time for each cycle

    protected:
        ifstream file;       // Source file object
        char *data_path;     // Path to source data
        char *fname_src;     // Name of source file
        char *data_type;     // Type of data in source ("scalar" or "vector")
        int nvals;           // 1 for scalar, 3 for vector
        int *dims;           // Dimensions of mesh for SILO output
        int Ntot;            // Product of dims elements
        int dims_in[ndims];  // Dimensions of source mesh (with HYM ghost zones)
        int Ntot_in;         // Product of dims_in elements
        long record_length;  // Binary record length for this data type
        char *stopmsg;       // Customizable error message
        char *varname;       // String name of variable (for SILO/VisIt)
        char *ascii_name;    // String prefix of output ASCII file (if written)

    public:
        HYMDataObj(char,char*,char*,int*,char*,int);
        ~HYMDataObj(void);
        static void ReadMesh_Binary(char*,char*,int*,float**,char*);
        virtual void WriteData_SILO(DBfile*,int,char*,float**) = 0;
        virtual void WriteData_ASCII(char*,int,double,float**) = 0;
        
    protected:
        static void StripGhostCoords(double*,float*&,int,int,int);
        static bool CompareMeshDims(int*,int*,char*,char*);
        void ValidateSourceFile(void);
        void ReadDims(long,int*);
        void GetTimes(void);
        void PositionPointer_Binary(int);
        void ReadVar_Binary(float*&);   
};

//============================================================================//
class HYMScalarObj : public HYMDataObj {
    public:
        HYMScalarObj(char,char*,char*,int*,char*);
        void WriteData_SILO(DBfile*,int,char*,float**);
        void WriteData_ASCII(char*,int,double,float**);
        
    protected:
        void ReadScalar_Binary(int,float*&);
};

//============================================================================//
class HYMVectorObj : public HYMDataObj {
    protected:
        char *varnames[ndims];

    public:
        HYMVectorObj(char,char*,char*,int*,char*);
        void WriteData_SILO(DBfile*,int,char*,float**);
        void WriteData_ASCII(char*,int,double,float**);
        
    protected:
        void ReadVector_Binary(int,float**);
};

//============================================================================//
//============================================================================//
