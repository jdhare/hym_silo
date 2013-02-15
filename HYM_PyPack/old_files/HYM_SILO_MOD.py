"""

Clayton Myers
HYM_SILO_MOD.py
Created:   9 July 2009
Modified: 11 December 2009

Module library of functions in support of the HYM_SILO.py conversion utility.

"""
#=============================================================================#
#=============================================================================#

import os
import sys
import time

#=============================================================================#
#=============================================================================#
def Execute_Operations(flags,paths,args):
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
    
    Print_Divider(level=1)
    
    #-------------------------------------------------------------------------#
    # Preparation tasks for the HYM-to-SILO conversion utility:
    if Convert_To_SILO:
        Locate_Data_Files(run_path,data_path,code_path)
        if not Locate_Directory(silo_path,"SILO output directory","  "):
            Get_User_Input("\n  Would you like to create the SILO output " + 
                           "directory?")
            Make_Directory(silo_path,"SILO output directory")
        Get_And_Print_Ncyc(data_path)
        Print_Divider(level=2)
    if Extract_Probe_Data:
        if not Locate_Directory(probe_path,"Probe output directory","  "):
            Make_Directory(probe_path,"Probe output directory")

    #-------------------------------------------------------------------------#
    # Start the timer:
    start_time = time.time()
    # Execute the C++ HYM-to-SILO conversion utility:
    if Convert_To_SILO:
        print "  Converting HYM Data to SILO Databases:"
        args = "%s %s %d %s" % (data_path,silo_path,cycle_arg,data_flags)
        os.system(r"./HYM_SILO.exe %s" % args)
        print "  Conversion operations completed."
        if Extract_Probe_Data:
            Print_Divider(level=2)
    # Execute the probe data extraction C++ utility:
    if Extract_Probe_Data:
        print "  Extracting the probe data:"
        os.system(r"./Probe_SILO.exe %s %s %d %f" % 
                    (silo_path,probe_path,cycle_arg,phi_rot))
        print "  Data extraction operations completed."
    # Output the timer results:
    Print_Divider(level=2)
    Print_Time(start_time,time.time())

    #-------------------------------------------------------------------------#
    Print_Divider(level=1)

#=============================================================================#
# def Extract_Probe_Data:
#    Make_Directory(probe_path)
#    Execute the C++ probe data extraction utility:
#    Print_Divider(level=2)
#    print "  Extracting Probe Data:"
#    outstr = r"./Probe_SILO.exe %s %s %d" % (silo_path,probe_path,cycle_arg)
#    os.system(outstr)
    
#=============================================================================#
#=============================================================================#
# Lists of raw HYM files required for conversion:
essential_files = ["hstat.d", "hgrid.d"]
nonessential_files = ["hmen.d", "history.d", "start.d", "0.dmp", "hfh3.d", 
                      "hfh3.l", "hfh3_ssx.d", "h3ds_ff.d", "t.err",  "t.out"]
del_files = ["0.dmp"]
copy_files = ["hybm.i","hmin.i","pscript"] # Fix this later (10/28/2009)

def filter_files(fname):
    """Function to filter the dirlist for raw HYM data files"""
    if fname in essential_files or fname in nonessential_files:
        return True
    else:
        return False

#=============================================================================#
def Locate_Data_Files(run_path,data_path,code_path):
    """
    Determines the location of the raw binary HYM data files for the run.  If 
    they are already stored in data_path, execution continues.  If they are not
    in data_path but new run data files exist in code_path, the user is 
    prompted as to whether these new files should be moved to data_path.
        
    """
    stored_files = []
    indent = "  "
    if Locate_Directory(run_path,"run archive directory",indent):
        if Locate_Directory(data_path,"data storage subdirectory",indent):
            dirlist = os.listdir(data_path)
            stored_files = filter(filter_files,dirlist)
            if len(stored_files) == 0:
                print(indent + "No stored files were found.")
            else:
                missing_files = []
                for item in essential_files:
                    if not item in stored_files:
                        missing_files.append(item)
                if len(missing_files) != 0:
                    print(indent + "An incomplete data set was located.")
                    Print_Divider(level=2)
                    print("  Files are missing from the stored data set:")
                    print "      Missing Files: ", missing_files
                    print("  Stopping script execution.")
                    Print_Divider(level=1)
                    sys.exit(0)
    if len(stored_files) == 0:
        Print_Divider(level=2)
        out_str = ("A suitable HYM data set was not found in storage.  Look " +
                   "in the active code directory for a new HYM data set?")
        Get_User_Input(out_str)
        print("")
        Move_Files(code_path,run_path,data_path)
        Print_Divider(level=2)

#=============================================================================#
def Move_Files(code_path,run_path,data_path):
    """
    Function to move the data files from the code directory init_path to their 
    storage destination in new_path.
        
    """
    if Locate_Directory(code_path,"active code directory","  "):
        # Assemble data file list to be moved
        dirlist = os.listdir(code_path)
        new_files = filter(filter_files,dirlist)
    
        # Move each file in the filtered dirlist to new_path from init_path
        move = len(new_files) != 0
        copy = False # "psi.dat" in dirlist
        delete = False # "0.dmp" in dirlist
        if move:
            # Ask whether the files should be moved:
            Print_Lines("  A new HYM data set was found in the active " +
                        "code directory.  By answering yes to this prompt, " +
                        "you are asking to do the following:","  ")
            print("      (1) Create the run archive directory " + 
                  "(if necessary).")
            print("      (2) Create the data storage subdirectory " + 
                  "(if necessary).")
            print("      (3) Move the new HYM data files to the data " +
                  "storage subdirectory.")
            if delete:
                print("      (4) Delete the file \"0.dmp\" if it exists.")
            Get_User_Input("Proceed with these actions?")
            
            # Make the storage directories if they do not exist:
            Make_Directory(run_path,"run archive directory")
            Make_Directory(data_path,"data storage subdirectory")
            # Employ overwrite protection for the new_path directory:
            Check_Overwrite(data_path,new_files)            
            # Move the files:
            print("\n  Moving data files to the directory \"%s\":" % data_path)
            for fname in new_files:
                print("    Moving file \"%s\" ..." % fname),
                fmove(fname,code_path,data_path)
                print("done.")
            if copy:
                fname = "psi.dat"
                print("    Copying file \"%s\" ..." % fname),
                fcopy(fname,code_path,data_path)
                print("done.")
            if delete:
                print("    Deleting file \"0.dmp\" ..."),
                Remove_File(os.path.join(code_path,"0.dmp"),message=False)
                print("done.")
            print("  Data file transfer complete.")
        else:
            print("  No new HYM data files were located.")
            Error_Exit("No suitable HYM data set was located.")        
    else:
        Error_Exit("The active code directory was not found.")


#=============================================================================#
def Check_Overwrite(path,flist):
    """
    Protects the data's destination directory from overwrite by requiring user 
    input in order to overwrite existing data files.
       
    """
    dirlist = os.listdir(path)
    conflict_list = []
    for fname in flist:
        if fname in dirlist:
            conflict_list.append(fname)
        
    if len(conflict_list) != 0:
        print("\n  File overwrite warning:" + 
              "\n      Destination directory: \"%s\"\n" % path)
        Get_User_Input("Should the files be overwritten?") 

#=============================================================================#
def fmove(fname,init_path,new_path):
    """UNIX move command"""
    init_pathname = os.path.join(init_path,fname)
    new_pathname = os.path.join(new_path,fname)
    os.system("mv %s %s" % (init_pathname,new_pathname))

def fcopy(fname,init_path,new_path):
    """UNIX copy command"""
    init_pathname = os.path.join(init_path,fname)
    new_pathname = os.path.join(new_path,fname)
    os.system("cp %s %s" % (init_pathname,new_pathname))

#=============================================================================#
def Locate_Directory(path,direc_str,indent):
    # print(indent + "Locating the %s." % direc_str)
    # print(indent + "Path: %s" % path)
    if os.path.exists(path):
        return True
    else:
        print(indent + "The %s was NOT FOUND." % direc_str)
        print(indent + "Path: %s" % path)
        return False

def Check_Directory(path):
    if not os.path.exists(path):
        out_str = "  The path \"%s\" does not exist.  \n" % path            + \
                  "  Stopping SILO file conversion."
        print(out_str)
        Print_Divider(level=1)
        sys.exit(0)    
    
def Make_Directory(path,direc_str):
    if not os.path.exists(path):
        os.mkdir(path)
        print("\n  The %s was created." % direc_str)
        print("  Path: %s" % path)

#=============================================================================#
def Open_File(fstring,permissions):
    """Open a file with error handling."""
    try:
        return open(fstring,permissions)
    except:
        print("\n  Unable to open the file \"%s\".\n" % fstring)
        sys.exit(0)      

def Remove_File(fstring,message=True):
    """Remove a file from the disk with error handling."""
    try: 
        os.remove(fstring)
        if(message):
            print("  Removing file \"%s\" from disk." % fstring)
    except: 
        print("\n  File \"%s\" not found for removal.\n" % fstring)

#=============================================================================#
def Get_User_Input(out_str):
    Print_Lines(out_str,"  ")
    print "      [y] Yes (continue)\n      [n] No  (exit)"
    in_str = raw_input("  Your selection is [n]: ")
    while in_str not in ["y","n",""]:
        in_str = raw_input("  Enter y or n as your selection [n]: ")
    if in_str == "y":
        return
    elif in_str == "n" or in_str == "":
        print("\n  Stopping script execution based on user input.")
        Print_Divider(level=1)
        sys.exit(0)
    else:
        raise ValueError

#=============================================================================#
def Get_And_Print_Ncyc(data_path):
    """Prints the number of cycles (Ncyc) for the run."""
    # Read words from status file 'hstat.d':
    fp = Open_File(os.path.join(data_path,"hstat.d"),"r")
    hstat_words = fp.read().split()
    fp.close()

    # Find number of cycles written during the run:
    i = hstat_words.index("i3dbout=")
    Ncycles = int(hstat_words[i+1])
    print("  Ncyc = %d" % Ncycles)

#=============================================================================#
def Print_Time(start_time,stop_time):
    t_tot = stop_time - start_time
    secs = t_tot%60
    mins = (t_tot-secs)/60
    print("  Total time elapsed: %0.2d:%05.2f" % (mins,secs))    

#=============================================================================#
def Print_Divider(level):
    """Prints the desired divider.  Level argument gives different divider."""
    if level == 1:
        print("\n=========================================================" + \
              "=======================\n")
    elif level == 2:
        print("\n  -------------------------------------------------------" + \
              "-----------------------\n")
    else:
        raise ValueError, level

#=============================================================================#
def Error_Exit(message):
    Print_Divider(level=2)
    print("  Error: %s" % message)
    print("  Stopping script execution.")
    Print_Divider(level=1)
    sys.exit(0)    

#=============================================================================#
def Print_Lines(string,indent,subindent=""):
    """Break up long strings into properly indented printable lines"""
    lines = []
    words = string.split(' ')
    while len(words) != 0:
        line = indent
        if len(lines) != 0:
            line += subindent
        lineinit = line
        while len(words) > 0 and len(line) + len(" "+words[0]) < 80:
            if line != lineinit:
                line += " "
            line += words[0]
            words.pop(0)
        if len(words) > 0 and len(words[0]) > 40:
            pass
        lines.append(line)
    for line in lines:
        print(line)    

#=============================================================================#
