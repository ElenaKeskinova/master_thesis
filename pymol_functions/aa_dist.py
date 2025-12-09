from pymol import cmd
import itertools

def residue_pairwise_distances(selection='my_selection', atom_name='CA', chain = "G"):
    # Get model (all atoms)
    model = cmd.get_model(selection)

    # Extract unique residues by CA atom
    residues = set()
    for atom in model.atom:
        # if atom.name == atom_name:
            res_id = (atom.chain, atom.resi, atom.resn)
            residues.add(res_id)
    
    for r in residues:
        print(f"{r[2]}{r[1]}")

    # Calculate distances between all pairs
    with open("distresult.txt", 'w') as f:
        f.write('Residue1,Residue2,Distance\n')
        
        for res1, res2 in itertools.combinations(residues, 2):
        
            dist_obj = f"temp_avg_{chain}_{res1[1]}_{chain}_{res2[1]}"
            d = cmd.distance(dist_obj,f"{chain} and resi {res1[1]}", f"{chain} and resi {res2[1]}")
            label1 = f"{res1[2]}{res1[1]}{res1[0]}"
            label2 = f"{res2[2]}{res2[1]}{res2[0]}"
            f.write(f'{label1},{label2},{d:.2f}\n')
            cmd.delete(dist_obj)

residue_pairwise_distances("antigen_contacts",chain = "antigen")