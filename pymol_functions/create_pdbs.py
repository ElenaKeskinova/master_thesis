from pymol import cmd
import os
import gc

file = "//BioInfo3/homes/Elena/pdb_structures/high_throughput_dock/b12/b12_controls.txt"
infolder = "//BioInfo3/homes/Elena/pdb_structures/high_throughput_dock/b12/peps/cntr"


with open(file, "r") as f:
    lines = f.readlines()
    n = len(lines)

for i, seq in enumerate(lines, start=1):
    seq = seq.strip()
    

    # Build peptide object
    obj_name = seq
    cmd.fab(seq, obj_name)

    # Save as PDB
    cmd.save(f"{infolder}/{obj_name}.pdb", obj_name)
    cmd.delete("all")
    if i % 100 == 0:  # every 50 cycles
        cmd.reinitialize()
    
print("ready")