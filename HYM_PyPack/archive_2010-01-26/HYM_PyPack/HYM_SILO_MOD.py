"""

Clayton Myers
HYM_SILO_MOD.py
Created:   9 July 2009
Modified: 11 December 2009

Module library of functions in support of the HYM_SILO.py conversion utility.

"""
#=============================================================================#
#=============================================================================#

import os, sys, time
import BASIC_FUNC_MOD as BFM
import HYM_FILES_MOD as HFM

#=============================================================================#
#=============================================================================#

essential_files = ["hstat.d", "hgrid.d"]

#=============================================================================#
def Convert_to_SILO(flags,paths,args):
    """Executes the utility operations based on user input"""
    Convert_To_SILO,Extract_Probe_Data = flags
    code_path,run_path,data_path,silo_path,probe_path = paths
    cycle_arg,data_flags,phi_rot = args
    
    if data_flags[0] == "1":
        essential_files.append("h3ds.d")
    if data_flags[1] == "1":
        essential_files.append("h3ds_ff.d")
    if data_flags[2] == "1":
        essential_files.append("h3db.d")
    if data_flags[3] == "1":
        essential_files.append("h3dv.d")
    if data_flags[4] == "1":
        essential_files.append("h3dj.d")
    
    BFM.Print_Divider(level=1)
    
    #-------------------------------------------------------------------------#
    # Preparation tasks for the HYM-to-SILO conversion utility:
    if Convert_To_SILO:
        Locate_Data_Set(run_path,data_path,code_path)
        if not BFM.Locate_Directory(silo_path,"SILO output directory","  "):
            BFM.Get_User_Input("\n  Would you like to create the SILO output " + 
                               "directory?")
            BFM.Make_Directory(silo_path,"SILO output directory")
            BFM.Print_Divider(level=2)
        Get_And_Print_Ncyc(data_path)
        BFM.Print_Divider(level=2)
    if Extract_Probe_Data:
        if not BFM.Locate_Directory(probe_path,"Probe output directory","  "):
            BFM.Make_Directory(probe_path,"Probe output directory")

    #-------------------------------------------------------------------------#
    # Start the timer:
    start_time = time.time()
    # Execute the C++ HYM-to-SILO conversion utility:
    if Convert_To_SILO:
        print "\n  Converting HYM Data to SILO Databases:"
        args = "%s %s %d %s" % (data_path,silo_path,cycle_arg,data_flags)
        os.system(r"./HYM_SILO.exe %s" % args)
        print "  Conversion operations completed."
        if Extract_Probe_Data:
            BFM.Print_Divider(level=2)
    # Execute the probe data extraction C++ utility:
    if Extract_Probe_Data:
        print "  Extracting the probe data:"
        os.system(r"./Probe_SILO.exe %s %s %d %f" % 
                    (silo_path,probe_path,cycle_arg,phi_rot))
        print "  Data extraction operations completed."
    # Output the timer results:
    BFM.Print_Divider(level=2)
    BFM.Print_Time(start_time,time.time())

    #-------------------------------------------------------------------------#
    BFM.Print_Divider(level=1,feedline=True)

#=============================================================================#
#=============================================================================#
def Locate_Data_Set(run_path, data_path, code_path):
    """
    Function to locate a HYM data set for SILO conversion. This function first
    checks in data_path.  If a data set is not located there but a new data set
    exists in code_path, the user is prompted as to whether these new files 
    should be moved to data_path.  The function also offers to create the run 
    archive directory at run_path and the data storage subdirectory at 
    data_path.
    
    """
    # First check data_path for an existing data set:
    if Validate_Data_Set(data_path):
        return
    # If not, ask to check in code_path:
    out_str = ("\n  A suitable HYM data set was not found in storage.  Look " +
               "in the active code directory for a new HYM data set?")
    BFM.Get_User_Input(out_str)
    # First locate the active code directory:
    if not BFM.Locate_Directory(code_path,"active code directory","  "):
        BFM.Error_Exit("The active code directory was not found.")
    # Check for a valid data set in code_path:
    if Validate_Data_Set(code_path):
        # Get user permission to make the file transfer:
        BFM.Print_Divider(level=2,feedline=True)
        BFM.Print_Lines("  A new HYM data set was found in the active " +
                        "code directory.  By answering yes to this " +
                        "prompt, you are asking to do the following:","  ")
        print("      (1) Create the run archive directory (if necessary).")
        print("      (2) Create the data storage subdirectory (if necessary).")
        print("      (3) Move the new HYM data files to the data storage " +
              "subdirectory.")
        BFM.Get_User_Input("Proceed with these actions?")       
        # Make the storage directories if they do not exist:
        BFM.Make_Directory(run_path,"run archive directory")
        BFM.Make_Directory(data_path,"data storage subdirectory")
        # Make the file transfer:
        HFM.File_Transfer(code_path,data_path,key_str="run_xfer",quiet=True)
    else:
        Error_Exit("No valid HYM data set could be located.")

#=============================================================================#
def Validate_Data_Set(path):
    """
    Function to determine if the directory at path holds a HYM data set that
    contains all essential files for HYM-to-SILO conversion.
        
    """
    indent = "  "
    # Locate and filter the data set:
    stored_files = []
    if BFM.Locate_Directory(path,"data directory",indent):
        dirlist = os.listdir(path)
        missing_files = []
        for fname in essential_files:
            if not fname in dirlist:
                missing_files.append(fname)
        if len(missing_files) == 0:
            return True
        elif len(missing_files) == len(essential_files):
            return False
        else:
            print(indent + "An incomplete data set was located.")
            BFM.Print_Divider(level=2)
            print("  Files are missing from the stored data set:")
            print "      Missing Files: ", missing_files
            BFM.Error_Exit("An incomplete data set was found.")
    else:
        return False

#=============================================================================#
def Get_And_Print_Ncyc(data_path):
    """Prints the number of cycles (Ncyc) for the run."""
    # Read words from status file 'hstat.d':
    fp = BFM.Open_File(os.path.join(data_path,"hstat.d"),"r")
    hstat_words = fp.read().split()
    fp.close()

    # Find number of cycles written during the run:
    if not "i3dbout=" in hstat_words:
        BFM.Error_Exit("Information is missing from the status file hstat.d.")
    i = hstat_words.index("i3dbout=")
    Ncycles = int(hstat_words[i+1])
    
    print("\n  Ncyc = %d" % Ncycles)

#=============================================================================#
#=============================================================================#
