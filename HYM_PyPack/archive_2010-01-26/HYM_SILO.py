"""

Clayton Myers
HYM_SILO.py
Created:  07 July 2009
Modified: 14 January 2010

This file is a wrapper script for the HYM-to-SILO conversion utility 
(HYM_SILO.exe) and the probe data extraction utility (Probe_SILO.exe).

User parameters:
    (1) Utility Flags (boolean):
          (a) Convert_Run_To_SILO -- Request HYM-to-SILO conversion.
          (b) Extract_Probe_Data  -- Request probe data extraction.
    (2) Paths (strings with leading and trailing slashes):
          (a) code_path  -- Path to the HYM code directory (used for moving 
                            data files from a new run to data_path if needed).
          (b) run_path   -- Path to the parent run storage directory.
    (3) cycle_arg (integer):
          -- An integer that determines which cycles will be converted by the 
             conversion utilities.  
          -- The term cycle refers to the sequential value of a given output 
             frame.  The number of cycles Ncyc corresponds to the integer value
             of "i3dbout" in the "hstat.d" run output file.
          -- If the cycle argument cycle_arg is between [1,Ncyc], then the 
             data will be converted for only that cycle value.
          -- If the cycle argument cycle_arg is zero (0), then all available
             cycles in the range [1,Ncyc] will be converted.
    (4) data_flags (string of 1's and 0's):
          -- A string of 1's and 0's that determines which variables will be
             written to the silo databases (e.g., 10110).
          -- The present data_flags string has five entries corresponding to:
                #1 = p (pressure)
                #2 = n (density)
                #3 = B (magnetic field)
                #4 = v (fluid velocity)
                #5 = J (current density)
    (5) phi_rot (float):
          -- Command line argument sent to Probe_SILO.exe that rotates the 
             synthesized probe array by phi_rot (in degrees).
"""

#==============================================================================#
#==============================================================================#
# User input variables:

Convert_To_SILO = True
Extract_Probe_Data = False

# Convert_To_SILO = False
# Extract_Probe_Data = True

# code_path = "./HYM_CoH/"
# run_path  = "./RunData_CoH/2009_07_21_co1/"
# code_path = "../scratch2/HYM_CounterH/"
# run_path  = "../scratch2/RunData_CounterH/2010_01_11_ct1/"

home_path        = "/global/homes/c/cmyers/Franklin/"
scratch_path     = home_path + "fscratch2/"
code_path = scratch_path + "HYM_CounterH/"
run_path  = scratch_path + "RunData_CounterH/2010_01_26_ct6/"

cycle_arg  = 0
data_flags = "10110"
phi_rot = 0.0

#==============================================================================#
################################################################################
################################################################################
#==============================================================================#
# Construct the derived paths and argument variables:
data_path  = run_path + "HYM_Data/" 
silo_path  = run_path + "SILO/"
probe_path = run_path + "Probe_Data/"

flags = [Convert_To_SILO,Extract_Probe_Data]
paths = [code_path,run_path,data_path,silo_path,probe_path]
args  = [cycle_arg,data_flags,phi_rot]

# Execute the script functionality:
import HYM_PyPack as HPP
HPP.Convert_to_SILO(flags,paths,args)

#==============================================================================#
################################################################################
################################################################################
#==============================================================================#
