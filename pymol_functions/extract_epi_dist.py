from pymol import cmd
import itertools

def analyze_contacts(pdb_path, antibody_chains='H+L', antigen_chains='G', cutoff=5.0, output_file='contact_distances.csv'):
    # Load the PDB
    cmd.reinitialize()
    cmd.load(pdb_path, 'complex')

    # Define selections
    cmd.select('antibody', f'chain {antibody_chains}')
    cmd.select('antigen', f'chain {antigen_chains}')

    # Find antigen atoms within cutoff of antibody
    cmd.select('antigen_contacts', f'antigen within {cutoff} of antibody')

    # Extract unique contact residues 
    model = cmd.get_model('antigen_contacts')
    contact_residues = []
    seen = set()

    for atom in model.atom:
        res_id = (atom.chain, atom.resi, atom.resn, atom.index)
        if (atom.chain, atom.resi) not in seen:
            seen.add((atom.chain, atom.resi))
            contact_residues.append(res_id)

    for r in contact_residues:
       print(f"{r[2]}{r[1]}")
    
    # Compute pairwise distances between atoms of contact residues
    with open(output_file, 'w') as f:
        f.write('Residue1,Residue2,Distance\n')
        for res1, res2 in itertools.combinations(contact_residues, 2):
            dist_obj = f"temp_avg_{res1[0]}_{res1[1]}_{res2[0]}_{res2[1]}"
            d = cmd.distance(dist_obj,f"chain {res1[0]} and resi {res1[1]}", f" chain {res2[0]} and resi {res2[1]}")
            label1 = f'{res1[2]}{res1[1]}{res1[0]}'
            label2 = f'{res2[2]}{res2[1]}{res2[0]}'
            f.write(f'{label1},{label2},{d:.2f}\n')
            cmd.delete(dist_obj)

    print(f"Pairwise distances written to: {output_file}")

# for pdb in pdbfiles:
#	analyze_contacts(pdb, output_file = f"{pdb}_cutoff5.txt", cutoff = 5)

def analyze_contacts2(antibody_chains='H+L', antigen_chains='G', cutoff=5.0, output_file='contact_distances.csv'): # for loaded model
    
    # Define selections
    cmd.select('antibody', f'chain {antibody_chains}')
    cmd.select('antigen', f'chain {antigen_chains}')

    # Find antigen atoms within cutoff of antibody
    cmd.select('antigen_contacts', f'antigen within {cutoff} of antibody')

    # Extract unique contact residues 
    model = cmd.get_model('antigen_contacts')
    contact_residues = []
    seen = set()

    for atom in model.atom:
        res_id = (atom.chain, atom.resi, atom.resn, atom.index)
        if (atom.chain, atom.resi) not in seen:
            seen.add((atom.chain, atom.resi))
            contact_residues.append(res_id)

    for r in contact_residues:
       print(f"{r[2]}{r[1]}")
    
    # Compute pairwise distances between atoms of contact residues
    with open(output_file, 'w') as f:
        f.write('Residue1,Residue2,Distance\n')
        for res1, res2 in itertools.combinations(contact_residues, 2):
            dist_obj = f"temp_avg_{res1[0]}_{res1[1]}_{res2[0]}_{res2[1]}"
            d = cmd.distance(dist_obj,f"chain {res1[0]} and resi {res1[1]}", f" chain {res2[0]} and resi {res2[1]}")
            label1 = f'{res1[2]}{res1[1]}{res1[0]}'
            label2 = f'{res2[2]}{res2[1]}{res2[0]}'
            f.write(f'{label1},{label2},{d:.2f}\n')
            cmd.delete(dist_obj)

    print(f"Pairwise distances written to: {output_file}")
