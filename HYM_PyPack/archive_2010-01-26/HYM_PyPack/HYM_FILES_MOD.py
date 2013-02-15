"""

Clayton Myers
HYM_FILES_MOD.py
Created: 13 January 2010

Module library of functions in support HYM file management.

"""
#==============================================================================#
#==============================================================================#

import os
import BASIC_FUNC_MOD as BFM

#==============================================================================#
#==============================================================================#
# Basic file lists:
data_ls  = ["hgrid.d", "h3ds.d", "h3ds_ff.d", "h3db.d", "h3dv.d", "h3dj.d"]
part_ls  = ["hstat.d", "hmen.d", "history.d", "hfh3.d", "hfh3_ssx.d", "hfh3.l"]
out_ls   = ["t.err", "t.out"]
cp_basic = ["hybm.p", "frcin_read.p", "hybm.i", "hmin.i", "fscript", "makefile"]

part_ls.extend(["h2ds_psi.d","h2ds_bt.d"])
data_ls.append("0.dmp")

# Other parameters:
indent = "      "

#==============================================================================#
#==============================================================================#
def File_Transfer(source_path, destination_path, key_str="manual", 
                  mv_ls=[], cp_ls=[], rm_ls=[], quiet=False):
    """Function to execute the desired file transfer"""    
    # Assign the proper values to the move/copy/remove lists:
    if key_str == "manual":
        xfer_str = "MANUAL"
    elif key_str == "run_xfer":
        mv_ls = data_ls[:] + part_ls[:] + out_ls[:]
        cp_ls = cp_basic[:]
        rm_ls = ["start.d"] # ,"0.dmp"]
        xfer_str = "RUN TRANSFER"
    elif key_str == "full_copy":
        cp_ls = data_ls[:] + part_ls[:] + out_ls[:] + cp_basic[:]
        xfer_str = "FULL COPY"
    elif key_str == "data_copy":
        cp_ls = data_ls[:] + part_ls[:]
        xfer_str = "DATA COPY"
    elif key_str == "part_copy":
        cp_ls = part_ls[:] + out_ls[:] + cp_basic[:]
        xfer_str = "PARTIAL COPY"
    elif key_str == "full_xfer":
        mv_ls = data_ls[:] + part_ls[:] + out_ls[:] + cp_basic[:]
        xfer_str = "FULL TRANSFER"
    elif key_str == "part_xfer":
        mv_ls = part_ls[:] + out_ls[:] + cp_basic[:]
        xfer_str = "PARTIAL TRANSFER"
    elif key_str == "purge":
        rm_ls = data_ls[:] + part_ls[:] + out_ls[:] + ["start.d"]
        xfer_str = "PURGE"
    else:
        BFM.Error_Exit("The key string \"%s\" is not recognized." % key_str)
    
    # Determine whether the destination directory is relevant:
    dest_needed = len(mv_ls) > 0 or len(cp_ls) > 0

    # Announce the file transfer:
    if not quiet:
        BFM.Print_Divider(level=1,feedline=True)
        print("  Initiating file transfer for HYM project files:")
        BFM.Print_Lines("Source Dir:      %s" % source_path,indent)
        if dest_needed:
            BFM.Print_Lines("Destination Dir: %s" % destination_path,indent)
        
    # Locate the source and directory:
    LD = BFM.Locate_Directory
    if not LD(source_path,"source directory",indent,quiet=True):
        BFM.Error_Exit("The source path for the file transfer does not exist.")

    # Locate the destination directory:
    dest_exists = LD(destination_path,"destination directory",indent,quiet=True)
    if dest_needed and not dest_exists:
        print("\n  The source directory has been located.")
        BFM.Get_User_Input("The destination directory was not found.  " +
                           "Create this directory?","  ")
        BFM.Make_Directory(destination_path,"destination directory","  ")
        
    # Match the file lists to the files in the source directory:
    dirlist = os.listdir(source_path)
    mv_filt = filter_files(mv_ls,dirlist)
    cp_filt = filter_files(cp_ls,dirlist)
    rm_filt = filter_files(rm_ls,dirlist)
    
    Nmv  = len(mv_ls)  ; Ncp  = len(cp_ls) ;  Nrm  = len(rm_ls)
    Nmvf = len(mv_filt); Ncpf = len(cp_filt); Nrmf = len(rm_filt)
 
    # Present the file transfer information:
    print("\n  Ready to execute the file transfer:")
    print(indent + "Transfer Description: %s" % xfer_str)
    print(indent + "Number of files to MOVE:   %02d/%02d" % (Nmvf,Nmv))
    print(indent + "Number of files to COPY:   %02d/%02d" % (Ncpf,Ncp))
    print(indent + "Number of files to REMOVE: %02d/%02d" % (Nrmf,Nrm))
    
    if (Nmvf+Ncpf+Nrmf) == 0:
        BFM.Error_Exit("No files were found to transfer.")
    else:
        if not quiet:
            BFM.Get_User_Input("Proceed with the file transfer?","  ") 
    
    # Check overwrite properties:
    BFM.Check_Overwrite(destination_path,mv_filt)
    BFM.Check_Overwrite(destination_path,cp_filt)
    
    # Execute the file transfer:
    print("\n  Executing the file transfer:")
    for fname in mv_filt:
        print("    Moving file \"%s\" ..." % fname),
        fmove(fname,source_path,destination_path)
        print("done.")
    for fname in cp_filt:
        print("    Copying file \"%s\" ..." % fname),
        fcopy(fname,source_path,destination_path)
        print("done.")
    for fname in rm_filt:
        print("    Deleting file \"%s\" ..." % fname),
        BFM.Remove_File(os.path.join(source_path,fname),message=False)
        print("done.")
    print("  Data file transfer complete.")    
    
    BFM.Print_Divider(level=1)

#==============================================================================#
def filter_files(itemlist,dirlist):
    """Function to filter the dirlist for raw HYM data files"""
    filterlist = []
    for fname in dirlist:
        if fname in itemlist:
            filterlist.append(fname)
    return filterlist
    
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

#==============================================================================#
#==============================================================================#
