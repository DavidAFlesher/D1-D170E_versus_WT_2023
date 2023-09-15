"""
This python script is meant to be loaded by ChimeraX through the command `runscript XXX` and
is intended to normalize PSII cryo-EM maps and perform important measurements. It requires
that the maps are already aligned within ChimeraX. The normalized voxel number is that which
gives the highest correlation coefficient  and those corresponding atom measurements should be 
used when doing distance comparisions. Voxel size is equivalent to pixel size.
"""


from chimerax.core.commands import run 
import numpy as np
import os, shutil, glob, time

def distance_measured(fiducial_model_int, fitting_model_int, fname):
    """
    To measure the distance between fiducial atoms within the EM map. These are ions which are typically high
    occupancy and thus allow more reliable centroid positions.

    Parameters
    ----------
    fiducial_model_int : int
        The model number within ChimeraX that contains the fiducial pdb.
        The atoms used are the non-heme iron, Mn1-4, OEC Ca, and Cl-1.
    fitting_model_int : int
        The model number within ChimeraX for the map which will be used for
        centroid position measurements.
    fname : str
        Name of the file to dump the measurements to. This is set with an f-string to contains
        the voxel number. Will contain distances between fiducial atoms.

    Returns
    -------
    None. Writes to file.

    """
    run(session, f"fitmap #{fiducial_model_int}@MN1 inMap #{fitting_model_int} moveWholeMolecules false gridStepMin 0.001")
    run(session, f"fitmap #{fiducial_model_int}@MN2 inMap #{fitting_model_int} moveWholeMolecules false gridStepMin 0.001")
    run(session, f"fitmap #{fiducial_model_int}@MN3 inMap #{fitting_model_int} moveWholeMolecules false gridStepMin 0.001")
    run(session, f"fitmap #{fiducial_model_int}@MN4 inMap #{fitting_model_int} moveWholeMolecules false gridStepMin 0.001")
    run(session, f"fitmap #{fiducial_model_int}@CA1 inMap #{fitting_model_int} moveWholeMolecules false gridStepMin 0.001")
    run(session, f"fitmap #{fiducial_model_int}@CL inMap #{fitting_model_int} moveWholeMolecules false gridStepMin 0.001")
    run(session, f"fitmap #{fiducial_model_int}@FE inMap #{fitting_model_int} moveWholeMolecules false gridStepMin 0.001")
    run(session, f"distance #{fiducial_model_int}@FE #{fiducial_model_int}@Cl ")
    run(session, f"distance #{fiducial_model_int}@FE #{fiducial_model_int}@MN1 ")
    run(session, f"distance #{fiducial_model_int}@FE #{fiducial_model_int}@MN2 ")
    run(session, f"distance #{fiducial_model_int}@FE #{fiducial_model_int}@MN3 ")
    run(session, f"distance #{fiducial_model_int}@FE #{fiducial_model_int}@MN4 ")
    run(session, f"distance #{fiducial_model_int}@FE #{fiducial_model_int}@CA1 ")
    run(session, f"distance #{fiducial_model_int}@MN1 #{fiducial_model_int}@MN2 ")
    run(session, f"distance #{fiducial_model_int}@MN1 #{fiducial_model_int}@MN3 ")
    run(session, f"distance #{fiducial_model_int}@MN1 #{fiducial_model_int}@MN4 ")
    run(session, f"distance #{fiducial_model_int}@MN1 #{fiducial_model_int}@CA1 ")
    run(session, f"distance #{fiducial_model_int}@MN2 #{fiducial_model_int}@MN3 ")
    run(session, f"distance #{fiducial_model_int}@MN2 #{fiducial_model_int}@MN4 ")
    run(session, f"distance #{fiducial_model_int}@MN2 #{fiducial_model_int}@CA1 ")
    run(session, f"distance #{fiducial_model_int}@MN3 #{fiducial_model_int}@MN4 ")
    run(session, f"distance #{fiducial_model_int}@MN3 #{fiducial_model_int}@CA1 ")
    run(session, f"distance #{fiducial_model_int}@MN4 #{fiducial_model_int}@CA1 ")
    run(session, f"save {fname}.pdb #{fitting_model_int}) # Optionally save the fitted pdb coordinates
    run(session, f"distance save {fname}.txt")
    run(session, "distance delete all")

def voxel_scale(input_voxel, model_int, ref_model_int, fiducial_model_int):
    """
    To set the voxel scale between two maps, with one being the reference to normalize
    against. Reminder, this requires the maps to be already nearly aligned and the map
    voxel size to be at high resolution, at most 1 A/pix. 

    Parameters
    ----------
    input_voxel : int
        The original voxel size of the map. Ex: 0.825 Å/pix
    model_int : int
        The model number within ChimeraX for the map which will be used for
        map to map fitting, or the map which will be compared against the reference map.
    ref_model_int : int
        The model number within ChimeraX for the reference map which will be used for
        map to map fitting used to evaluate the voxel scaling mismatch.
    fiducial_model_int : int
        The model number within ChimeraX that contains the fiducial pdb.
        The atoms used are the non-heme iron, Mn1-4, OEC Ca, and Cl-1.

    Returns
    -------
    None. Write to file.

    """
    if os.path.exists("logs"):
        shutil.rmtree("logs")
    time.sleep(3) # fix weird bug where ChimeraX wouldn't delete the dir before mkdir
    os.mkdir("logs")
    #pos_voxel_arr = np.arange(input_voxel, input_voxel+0.02, 0.001)
    neg_voxel_arr = np.arange(input_voxel-0.015, input_voxel+0.001, 0.001)
    full_voxel_arr = np.arange(input_voxel-0.015, input_voxel+0.015, 0.001)
    


    #prime the model. You need to slowly adjust to the bottom
    for voxel in neg_voxel_arr[::-1]: # reverse for high to low
        voxel = round(voxel, 3)
        run(session, f"fitmap #{str(model_int)} inMap #{str(ref_model_int)} metric overlap gridStepMin 0.01")
        run(session, f"volume #{str(model_int)} voxelSize {str(voxel)}")
        run(session, f"fitmap #{str(model_int)} inMap #{str(ref_model_int)} metric overlap gridStepMin 0.01")
        run(session, "wait 5")

    # now test the voxel scales
    for voxel in full_voxel_arr:
        voxel = round(voxel, 3)
        fname_fit_logs = os.path.join("logs", f"fit_logs_pixelSize_{voxel}")
        fname_distance_logs = os.path.join("logs", f"fit_distance_pixelSize_{voxel}")
        run(session, f"fitmap #{str(model_int)} inMap #{str(ref_model_int)} metric overlap gridStepMin 0.01 logFits {fname_fit_logs}.txt")
        run(session, f"volume #{str(model_int)} voxelSize {str(voxel)}")
        run(session, f"fitmap #{str(model_int)} inMap #{str(ref_model_int)} metric overlap gridStepMin 0.01 logFits {fname_fit_logs}.txt")
        run(session, "wait 5")
        # now measure my distances
        distance_measured(fiducial_model_int, model_int, fname_distance_logs)

def save_data2csv(csv_fname):
    # %%
    """
    Parse the data logs to an easily readable CSV file.

    Parameters
    ----------
    fiducial_model_int : str
        File name to save the CSV file as. Should contain the normalized map and the reference map names.
    Returns
    -------
    None. Save to file.
    """
    # Get log file names. 
    distance_files = glob.glob(os.path.join("logs", "fit_distance_*"))
    fit_files = glob.glob(os.path.join("logs", "fit_logs_*"))

    # %% Parse correlations

    # define dictionary that contains the voxel size as the key, then key:value for
    # different parameters, e.g. Fe - Mn1 distance, correlation, etc.
    data_dict = dict()
    for file in fit_files: 
        with open(file, 'r') as f:
            pixelSize = os.path.splitext(file)[0].split("_")[-1] # parse from fname
            lines = f.readlines()
            # the data is stored as plain text. parse the relevant info as follows.
            labels = lines[0].split()
            values = lines[-1].split()
            # add to dict as follows
            data_dict[pixelSize] = dict()
            for l, v in zip(labels, values):
                if l == "correlation":
                    data_dict[pixelSize][l] = v
                elif l == "correlation_about_mean":
                    data_dict[pixelSize][l] = v
                elif l == "overlap":
                    data_dict[pixelSize][l] = v

    # %% Parse distances
    for file in distance_files:
        with open(file, 'r') as f:
            pixelSize = os.path.splitext(file)[0].split("_")[-1] # parse from fname
            lines = f.readlines()
            # the data is stored as plain text. parse the relevant info as follows.
            for line in lines[3:]:
                line_split = line.split()
                atom1 = line_split[3]
                atom2 = line_split[-2].replace(":", "") # drop :
                distance = line_split[-1]
                distance_label = f"{atom1}-{atom2}_distance"
                data_dict[pixelSize][distance_label] = distance
                
    # %% Save data to csv
    voxels = list(data_dict.keys())
    data_headers = list(data_dict[voxels[-1]].keys())
    csv_header = f"voxel size,{','.join(data_headers)}\n"
    with open(f"{csv_fname}.csv", 'w') as f:
        f.write(csv_header)
        for voxel in voxels: # input from the data_dict corresponding to the header
            v = list(data_dict[voxel].values())
            write_line = f"{voxel},{','.join(v)}\n"
            f.write(write_line)


input_voxel = 0.825
model_int = 14
ref_model_int = 1
fiducial_model_int = 6
csv_fname = "voxels_D170E-jimin_to_D170E-jimin"

voxel_scale(input_voxel, model_int, ref_model_int, fiducial_model_int)
save_data2csv(csv_fname)
