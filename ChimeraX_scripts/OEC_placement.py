"""
This python script is meant to be loaded by ChimeraX through the command `runscript XXX` and
is intended to fit the OEC (and other) metal atoms at their signal maximum, then measure OEC distances.
"""


from chimerax.core.commands import run 

def distance_measured(pdb_model, map_model, fname):
    run(session, f"fitmap #{pdb_model}/A@MN1 inMap #{map_model} moveWholeMolecules false gridStepMin 0.001")
    run(session, f"fitmap #{pdb_model}/A@MN2 inMap #{map_model} moveWholeMolecules false gridStepMin 0.001")
    run(session, f"fitmap #{pdb_model}/A@MN3 inMap #{map_model} moveWholeMolecules false gridStepMin 0.001")
    run(session, f"fitmap #{pdb_model}/A@MN4 inMap #{map_model} moveWholeMolecules false gridStepMin 0.001")
    run(session, f"fitmap #{pdb_model}/A@CA1 inMap #{map_model} moveWholeMolecules false gridStepMin 0.001")
    run(session, f"fitmap #{pdb_model}/A:403@CL inMap #{map_model} moveWholeMolecules false gridStepMin 0.001")
    run(session, f"fitmap #{pdb_model}/A:404@CL inMap #{map_model} moveWholeMolecules false gridStepMin 0.001")
    run(session, f"fitmap #{pdb_model}/A:402@FE inMap #{map_model} moveWholeMolecules false gridStepMin 0.001")
    run(session, f"distance #{pdb_model}/A@MN1 #{pdb_model}/A@MN2 ")
    run(session, f"distance #{pdb_model}/A@MN1 #{pdb_model}/A@MN3 ")
    run(session, f"distance #{pdb_model}/A@MN1 #{pdb_model}/A@MN4 ")
    run(session, f"distance #{pdb_model}/A@MN1 #{pdb_model}/A@CA1 ")
    run(session, f"distance #{pdb_model}/A@MN2 #{pdb_model}/A@MN3 ")
    run(session, f"distance #{pdb_model}/A@MN2 #{pdb_model}/A@MN4 ")
    run(session, f"distance #{pdb_model}/A@MN2 #{pdb_model}/A@CA1 ")
    run(session, f"distance #{pdb_model}/A@MN3 #{pdb_model}/A@MN4 ")
    run(session, f"distance #{pdb_model}/A@MN3 #{pdb_model}/A@CA1 ")
    run(session, f"distance #{pdb_model}/A@MN4 #{pdb_model}/A@CA1 ")
    run(session, f"save {fname}.pdb #{pdb_model}") # Optionally save the fitted pdb coordinates
    run(session, f"distance save {fname}.txt")
    run(session, "distance delete all")
   

pdb_model = 2 # Model number assigned to the data in ChimeraX
map_model = 1 # Model number assigned to the data in ChimeraX
fname_distance_logs = "fit1_distances"
distance_measured(pdb_model, map_model, fname_distance_logs)

