//============================================================================//
/*

Clayton Myers
SILO_CycleObj.cpp
Created:  06 February 2010
Modified: 06 February 2010

Header file for SILO cycle objects.

*/
//============================================================================//
//============================================================================//
class SILO_CycObj {
    public:
        int cycle;               // Cycle number for this object
        double time;             // Time for this cycle
        char *stat_str;          // String for report status
        char *time_str;          // String for report time value
        bool report_flag;        // Flagged if this database requires a report
        bool write_flag;         // Flagged if this database should be written
        bool mask_flags[nvars];  // Mask array for writing each data member
        
    protected:
        int *silodims;               // Dimensions of mesh for SILO output
        float **mesh_coords;     // Coordinates of the HYM mesh
        char *stopmsg;           // Customizable error message
        HYMDataObj **data_objs;  // Vector of HYM data objects

    public:
        SILO_CycObj(int,float**,int*,HYMDataObj**,bool*,char*);
        ~SILO_CycObj(void);
        void Write_SILO(char*);
        void Write_ASCII(char*);
        
    protected:
        static void StripGhostCoords(double*,float*&,int,int,int);
        static bool CompareMeshDims(int*,int*,char*,char*);
        void SetFlags(bool*);
        void SetTime(void);
};

//============================================================================//
//============================================================================//
