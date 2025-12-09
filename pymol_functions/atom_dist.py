from pymol import cmd

# distances for each pair of atoms
def pairwise_dist(selection='all', cutoff=3.0):
    model = cmd.get_model(selection)
    atoms = model.atom
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            d = cmd.get_distance(f"id {atoms[i].index}", f"id {atoms[j].index}")
            if d <= cutoff:
                print(f"{atoms[i].resn}{atoms[i].resi}.{atoms[i].name} - {atoms[j].resn}{atoms[j].resi}.{atoms[j].name} : {d:.2f} Å")

pairwise_dist("antigen_contacts", 5.0)