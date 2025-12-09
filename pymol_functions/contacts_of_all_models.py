from pymol import cmd

def write_contacts(path= "1bj1_2_docked"):
    files = os.listdir(path)
    pdbfiles = [f for f in files if f.endswith(".pdb")]

    for pf in pdbfiles:
        # Load the PDB
        cmd.reinitialize()
        cmd.load(f"{path}/{pf}", 'complex')
    
        get_contacts(antibody='chain A', antigen='chain B', cutoff=5.0, output_file=f'{path}/{pf}_contacts.txt')
        print(f"{pf} ready")


