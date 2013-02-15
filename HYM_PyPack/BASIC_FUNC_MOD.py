"""

Clayton Myers
BASIC_FUNC_MOD.py
Created: 13 January 2010

Module library of basic functions in support HYM-to-SILO scripting operations.

"""
#==============================================================================#
#==============================================================================#

import os, sys, time

#==============================================================================#
#==============================================================================#
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
              "\n      Destination directory: %s" % path)
        Get_User_Input("Should the files be overwritten?") 

#==============================================================================#
def Locate_Directory(path,direc_str,indent,quiet=False):
    # print(indent + "Locating the %s." % direc_str)
    # print(indent + "Path: %s" % path)
    if os.path.exists(path):
        return True
    else:
        if not quiet:
            print ""
            Print_Lines("The %s was NOT FOUND." % direc_str,indent)
            Print_Lines("Path: %s" % path, indent)
        return False

def Check_Directory(path):
    if not os.path.exists(path):
        out_str = "  The path \"%s\" does not exist.  \n" % path            + \
                  "  Stopping SILO file conversion."
        print(out_str)
        Print_Divider(level=1)
        sys.exit(0)    
    
def Make_Directory(path,direc_str,indent="  "):
    if not os.path.exists(path):
        os.mkdir(path)
        print ""
        Print_Lines("The %s was created." % direc_str,indent)
        Print_Lines("Path: %s" % path,indent)

#==============================================================================#
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
        print("  File \"%s\" not found for removal." % fstring)
        
#==============================================================================#
def Get_User_Input(out_str,indent="  "):
    Print_Lines(out_str,indent)
    print indent + "    [y] Yes (continue)\n" + indent + "    [n] No  (exit)"
    in_str = raw_input(indent + "Your selection is [n]: ")
    while in_str not in ["y","n",""]:
        in_str = raw_input(indent + "Enter y or n as your selection [n]: ")
    if in_str == "y":
        return
    elif in_str == "n" or in_str == "":
        print("\n  Stopping script execution based on user input.")
        Print_Divider(level=1,feedline=True)
        sys.exit(0)
    else:
        raise ValueError

#==============================================================================#
def Print_Divider(level,feedline=False):
    """Prints the desired divider.  Level argument gives different divider."""
    if level == 1:
        print("\n=========================================================" + \
              "=======================")
    elif level == 2:
        print("\n  -------------------------------------------------------" + \
              "-----------------------")
    else:
        raise ValueError, level
    if feedline:
        print("")

#==============================================================================#
def Error_Exit(message):
    Print_Divider(level=2,feedline=True)
    Print_Lines("Error: %s" % message,"  ")
    print("  Stopping script execution.")
    Print_Divider(level=1,feedline=True)
    sys.exit(0)    

#==============================================================================#
def Print_Lines(string,indent,subindent=""):
    """Break up long strings into properly indented printable lines"""
    lines = []
    words = string.split(' ')
    maxiters = 7
    iters = 0
    while len(words) != 0:
        if iters > maxiters:
            Error_Exit("Infinite while loop hit in Print_Lines.")
        line = indent
        if len(lines) != 0:
            line += subindent
        lineinit = line
        while len(words) > 0 and (len(line) + len(" "+words[0])) < 80:
            if line != lineinit:
                line += " "
            line += words[0]
            words.pop(0)
        if len(line) < 80 and len(words) > 0 and len(words[0]) > 60:
            line += " "
            subindent += "    "
            while len(line) < 80:
                line += words[0][0]
                words[0] = words[0][1:]
        lines.append(line)
        iters += 1
    for line in lines:
        print(line)
        
#=============================================================================#
def Print_Time(start_time,stop_time):
    t_tot = stop_time - start_time
    secs = t_tot%60
    mins = (t_tot-secs)/60
    print("\n  Total time elapsed: %0.2d:%05.2f" % (mins,secs))      

#==============================================================================#
#==============================================================================#
