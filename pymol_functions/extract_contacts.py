from pymol import cmd


def get_contacts(antibody='antibody', antigen='ag', cutoff=5.0, output_file='paratope_contacts.csv'):
    # Load the PDB
    # cmd.reinitialize()
    # cmd.load(pdb_path, 'complex')

    # Define selections
    # cmd.select('antibody', f'chain {antibody_chains}')
    # cmd.select('antigen', f'chain {antigen_chains}')

    # Find antigen atoms within cutoff of antibody
    cmd.select('antigen_contacts', f'{antigen} within {cutoff} of {antibody}')

    # Extract unique contact residues 
    model = cmd.get_model('antigen_contacts')
    contact_residues = []
    seen = set()

    for atom in model.atom:
        res_id = (atom.chain, atom.resi, atom.resn, atom.index)
        if (atom.chain, atom.resi) not in seen:
            seen.add((atom.chain, atom.resi))
            contact_residues.append(res_id)

    
    # Compute pairwise distances between atoms of contact residues
    with open(output_file, 'w') as f:
        f.write('AG_contact\tAB_contacts\n')
        for ag_res in contact_residues:
            
            label1 = f'{ag_res[2]}{ag_res[1]}{ag_res[0]}'
            
            # for ab_contacts
            ab_list = []
            seen = set()
            
            cmd.select('ab_contacts', f'{antibody} within {cutoff} of chain {ag_res[0]} and resi {ag_res[1]}')
            ab_res = cmd.get_model('ab_contacts')
            for atom in ab_res.atom:
                res_id = (atom.chain, atom.resi, atom.resn, atom.index)
                if (atom.chain, atom.resi) not in seen:
                    seen.add((atom.chain, atom.resi))
                    ab_list.append(res_id)
            # print(ab_list)
            print(f"{len(ab_list)} contacts of {label1}")
            
            label2 = ",".join([f'{res[2]}{res[1]}{res[0]}' for res in sorted(ab_list)])
            
            f.write(f'{label1}\t{label2}\n')

    print(f"contact residues written to: {output_file}")