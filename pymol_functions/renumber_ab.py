from pymol import cmd

def renumber_antibody(pdb_path, heavy_chain='H', light_chain='L', output_pdb='renumbered_antibody.pdb'):
    # Start fresh and load structure
    cmd.reinitialize()
    cmd.load(pdb_path, "antibody")

    # Get number of residues in heavy chain
    heavy_residues = cmd.count_atoms(f"chain {heavy_chain} and name CA")
    if heavy_residues == 0:
        print(f"Error: No residues found in chain {heavy_chain}")
        return
    print(f"Heavy chain residues: {heavy_residues}")

    # Renumber heavy chain starting from 1
    cmd.alter(f"chain {heavy_chain}", "resi=str(int(resi)-int(resi)+1+int(resi)-int(resi))")
    cmd.alter(f"chain {heavy_chain}", "resi=str(int(resi)+enumerate(range(1,heavy_residues+1)))")

    # Get number of residues in light chain
    light_residues = cmd.count_atoms(f"chain {light_chain} and name CA")
    if light_residues == 0:
        print(f"Error: No residues found in chain {light_chain}")
        return
    print(f"Light chain residues: {light_residues}")

    # Renumber light chain to continue from heavy chain
    offset = heavy_residues
    counter = 1
    model = cmd.get_model(f"chain {light_chain} and name CA")
    for atom in model.atom:
        cmd.alter(f"chain {light_chain} and resi {atom.resi}", f"resi='{offset+counter}'")
        counter += 1

    # Update and save
    cmd.sort()
    cmd.save(output_pdb)
    print(f"Renumbered structure saved to {output_pdb}")