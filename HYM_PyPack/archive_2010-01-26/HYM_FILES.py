"""

Clayton Myers
HYM_FILES.py
Created: 13 January 2010

This script calls functionality from HYM_PyPack to move around files from a
HYM data set.  The user must set a source path, a destination path, and the
key string (key_str) that tells the file transfer function which type of file
transfer should be executed.

The following built-in key strings are available to the user:
    (1) run_xfer:  move all data files, copy source files, remove start.d, 0.dmp
    (2) full_copy: copy every file in the data set
    (3) part_copy: copy every file in the data set except high-storage *.d files
    (4) full_xfer: move every file in the data set
    (5) part_xfer: move every file in the data set except high-storage *.d files
    (6) manual:    user must supply the mv_ls, cp_ls, rm_ls lists of filename 
                   strings as arguments to the File_Transfer function.
"""

#==============================================================================#
#==============================================================================#
# User input variables:

key_str = "part_copy"
# key_str = "data_copy"
# key_str = "purge"

runID   = "2010_01_26_ct6"

#------------------------------------------------------------------------------#

home_path        = "/global/homes/c/cmyers/Franklin/"
scratch_path     = home_path + "fscratch2/"
code_path        = scratch_path + "HYM_CounterH/"
run_path_home    = home_path    + "RunData_CounterH/" + runID + "/"
run_path_scratch = scratch_path + "RunData_CounterH/" + runID + "/"

if key_str == "part_copy":
    source_path      = run_path_scratch + "HYM_Data/"
    destination_path = run_path_home    + "HYM_Data/"
elif key_str == "data_copy":
    source_path      = run_path_scratch + "HYM_Data/"
    destination_path = code_path
elif key_str == "purge":
    source_path      = code_path
    destination_path = ""
else:
    print "\n    Invalid key string.  Exiting program.\n"
    sys.exit(0)

#==============================================================================#
#==============================================================================#
# Execute the script functionality:

import HYM_PyPack as HPP
HPP.File_Transfer(source_path,destination_path,key_str)
print ""

#==============================================================================#
#==============================================================================#
