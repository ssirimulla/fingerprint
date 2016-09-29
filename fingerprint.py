###############
#
#  Based on PyPlif
#
###############

import getopt, sys, os, copy, glob
import numpy as np
from openbabel import OBMol, OBConversion, OBMolAtomIter, OBResidueIter, OBResidueAtomIter, OBAtomAtomIter
from optparse import OptionParser
from bitarray import bitarray
from collections import OrderedDict
from matplotlib import mpl,pyplot
from matplotlib.ticker import MultipleLocator, FixedLocator

protein_ref = ""
ligand_ref = ""
output = "fingerprint.txt"
graph_name = "fingerprint.tiff"
graph_title = ""

parser = OptionParser()
#Protein file
parser.add_option("-p", "--protein", dest="protein",
                  help="read config from FILE", metavar="FILE")
#Ligand file
parser.add_option("-l", "--ligand", dest="ligand",
                  help="read config from FILE", metavar="FILE")
#Output file
parser.add_option("-o", "--output", dest="output",
                  help="read config from FILE", metavar="FILE")
#Graph file name
parser.add_option("-g", "--graph_name", dest="graph_name",
                  help="read config from FILE", metavar="FILE")
#Title for graph
parser.add_option("-t", "--graph_title", dest="graph_title",
                  help="read config from FILE", metavar="FILE")

(options_parser, args) = parser.parse_args()

if options_parser.protein:
    protein_ref = options_parser.protein
if options_parser.ligand:
    ligand_ref = options_parser.ligand
if options_parser.output:
    output = options_parser.output
if options_parser.graph_name:
    graph_name = options_parser.graph_name
if options_parser.graph_title:
    graph_title = options_parser.graph_title

conv = OBConversion()
conv.SetInFormat("pdbqt")

protein  = OBMol()
ligand  = OBMol()

conv.ReadFile(protein, protein_ref)
conv.ReadFile(ligand, ligand_ref)

ligand_file = open(ligand_ref, 'r')
ligand_lines = [line for line in ligand_file]
ligand_atoms = []
ligand_charges = []
ligand_coords = []
ligand_box = [] #Box for each model that limits the possibility of interactions
i = 0
first_model = True
conformations = 0

for line in ligand_lines:
    newLine = ' '.join(line.split()).split(' ')
    if newLine[0] == "ENDMDL":
        xmax = max(ligand_coords, key=lambda x:x[0])[0]
        ymax = max(ligand_coords, key=lambda x:x[1])[1]
        zmax = max(ligand_coords, key=lambda x:x[2])[2]
        xmin = min(ligand_coords, key=lambda x:x[0])[0]
        ymin = min(ligand_coords, key=lambda x:x[1])[1]
        zmin = min(ligand_coords, key=lambda x:x[2])[2]
        i += 1
        conformations += 1
        first_model = False
        ligand_box.append((xmax,xmin,ymax,ymin,zmax,zmin))
        ligand_coords = []

    elif newLine[0] == "HETATM" and first_model:
        ligand_atoms.append(newLine[-1])
        ligand_charges.append(newLine[-2])
    if newLine[0] == "HETATM":
        ligand_coords.append((float(newLine[4]), float(newLine[5]), float(newLine[6])))




protein_file = open(protein_ref, 'r')
protein_lines = [line for line in protein_file]
resnum = ""
resname = []
holder = []
protein_atoms = []
protein_charges = []
res_atoms = []
res_charges = []
protein_coords = []
residue_dict = {} #stores bitarrays (the fingerprint output)
resid = []
i = -1
for line in protein_lines:
    newLine = ' '.join(line.split()).split(' ')
    if newLine[0] == "ATOM" and not int(line[23:26]) == resnum:
        protein_atoms.append([newLine[-1]])
        protein_charges.append([newLine[-2]])
        if i > 0:
            xmax = max(protein_coords, key=lambda x:x[0])[0]
            ymax = max(protein_coords, key=lambda x:x[1])[1]
            zmax = max(protein_coords, key=lambda x:x[2])[2]
            xmin = min(protein_coords, key=lambda x:x[0])[0]
            ymin = min(protein_coords, key=lambda x:x[1])[1]
            zmin = min(protein_coords, key=lambda x:x[2])[2]
            for k in range(len(ligand_box)):
                if (ligand_box[k][1]-4.5 <= xmax <= ligand_box[k][0]+4.5 or ligand_box[k][1]-4.5 <= xmin <= ligand_box[k][0]+4.5) and (ligand_box[k][3]-4.5 <= ymax <= ligand_box[k][2]+4.5 or ligand_box[k][3]-4.5 <= ymin <= ligand_box[k][2]+4.5) and (ligand_box[k][5]-4.5 <= zmax <= ligand_box[k][4]+4.5 or ligand_box[k][5]-4.5 <= zmin <= ligand_box[k][4]+4.5):
                    resname.append(resid[i])
                    res_atoms.append(protein_atoms[i])
                    res_charges.append(protein_charges[i])
                    break
            i += 1
            protein_coords = []
        else:
            i += 1
        resnum = int(line[23:26])
        resid.append(line[17:20]+str(resnum))
    elif newLine[0] == "ATOM":
        protein_atoms[i].append(newLine[-1])
        protein_charges[i].append(newLine[-2])
        protein_coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))

protein_rings = protein.GetSSSR()
pro_ring_dict = {}
pro_ring_num = 0
for ring in protein_rings:
    if ring.IsAromatic():
        pro_ring_dict[pro_ring_num] = [protein.GetAtom(ring._path[0]).GetResidue().GetName()+str(protein.GetAtom(ring._path[0]).GetResidue().GetNum()), ring._path]
        pro_ring_num += 1

print 'done reaing molecules, starting interaction analysis'
#Cross product of ring (to get the normal)
def getringcross(coords):
    a = coords[0]-coords[1]
    b = coords[0]-coords[2]
    crossprod = np.cross(a,b)
    return crossprod
#Checks if an atom in one ring is within 4A of an atom in another ring
def ringdistance(ring1, ring2, mol1, mol2):
    for atom1 in ring1:
        for atom2 in ring2._path:
            atomdistance = mol1.GetAtom(atom1).GetDistance(mol2.GetAtom(atom2))
            if atomdistance <= 4.0:
                return True
    return False
#Nonpolar based on SMARTS patterns
def isnonpolar(atom):
    if atom.MatchesSMARTS('[#6,#16,F,Cl,Br,I]'):
        return True
    return False


RADTODEG = 180 / np.pi
output_holder = []
for conform in range(conformations):
    if conform > 0:
        conv.Read(ligand)
    for residue in resname:
        residue_dict[residue] = bitarray('00000000')
    #Finda ring interactions
    for ring in pro_ring_dict:
        if pro_ring_dict[ring][0] in resname:
            pro_coords = []
            numcoord = 0
            while numcoord < 3:
                pro_coords.append(np.array([protein.GetAtom(pro_ring_dict[ring][1][numcoord]).x(),
                protein.GetAtom(pro_ring_dict[ring][1][numcoord]).y(), protein.GetAtom(pro_ring_dict[ring][1][numcoord]).z()]))
                numcoord += 1
            for lig_ring in ligand.GetSSSR():
                if lig_ring.IsAromatic():
                    lig_coords = []
                    numcoord = 0
                    while numcoord < 3:
                        lig_coords.append(np.array([ligand.GetAtom(lig_ring._path[numcoord]).x(),
                                   ligand.GetAtom(lig_ring._path[numcoord]).y(), ligand.GetAtom(lig_ring._path[numcoord]).z()]))
                        numcoord += 1
                    pro_ring = pro_ring_dict[ring][1]
                    in_range = ringdistance(pro_ring, lig_ring, protein, ligand)
                    if in_range:
                        pro_cross = getringcross(pro_coords)
                        lig_cross = getringcross(lig_coords)
                        dot = np.dot(pro_cross, lig_cross)
                        pro_modulus = np.sqrt((pro_cross*pro_cross).sum())
                        lig_modulus = np.sqrt((lig_cross*lig_cross).sum())
                        cos_angle = dot / pro_modulus / lig_modulus
                        ring_angle = np.arccos(cos_angle) * RADTODEG

                        if (30.0 >= ring_angle) or (150.0 <= ring_angle):
                            #Face to face
                            residue_dict[pro_ring_dict[ring][0]] |= bitarray('01000000')
                        if (30.0 <= ring_angle <= 150.0):
                            #Edge to face
                            residue_dict[pro_ring_dict[ring][0]] |= bitarray('00100000')

#    print 'done with ring interactions with conformation: '+str(conform+1)

    i = 0
    for residue in OBResidueIter(protein):
        residuename = residue.GetName() + str(residue.GetNum())
        if residuename in resname:
            i = resname.index(residuename)
            j = 0
            for res_atom in OBResidueAtomIter(residue):
                if not res_atom.IsHydrogen():
                    k = 0
                    for lig_atom in OBMolAtomIter(ligand):
                        if not lig_atom.IsHydrogen():
                            dist = res_atom.GetDistance(lig_atom)
                            if dist <= 4.5:
                                if isnonpolar(res_atom) and isnonpolar(lig_atom):
                                    #Apolar
                                    residue_dict[residuename] |= bitarray('10000000')
                                if dist <= 4.0:
                                    if float(res_charges[i][j]) > 0. and float(ligand_charges[k]) < 0.:
                                        #Electrostatic (protein positive)
                                        residue_dict[residuename] |= bitarray('00000100')
                                    if float(res_charges[i][j]) < 0. and float(ligand_charges[k]) > 0.:
                                        #Electrostatic (protein negative)
                                        residue_dict[residuename] |= bitarray('00000010')
                                    if dist <= 3.5:
                                        if res_atom.IsHbondDonor() and lig_atom.IsHbondAcceptor():
                                            #H-bond (protein donor)
                                            for neighborDon in OBAtomAtomIter(res_atom):
                                                if neighborDon.IsHydrogen():
                                                    angle = res_atom.GetAngle(neighborDon, lig_atom)
                                                    if angle > 135.0:
                                                        residue_dict[residuename] |= bitarray('00010000')
                                        if res_atom.IsHbondAcceptor() and lig_atom.IsHbondDonor():
                                            #H-bond (protein acceptor)
                                            for neighborDon in OBAtomAtomIter(lig_atom):
                                                if neighborDon.IsHydrogen():
                                                    angle = lig_atom.GetAngle(neighborDon, lig_atom)
                                                    if angle > 135.0:
                                                        residue_dict[residuename] |= bitarray('00001000')
                            #Halogens
                            if lig_atom.MatchesSMARTS('[Cl]') and res_atom.IsHbondAcceptor() and dist <= 3.5:
                                for neighborHal in OBAtomAtomIter(lig_atom):
                                    if neighborHal.IsCarbon():
                                        angle = neighborHal.GetAngle(lig_atom, res_atom)
                                        if angle > 146.0:
                                            residue_dict[residuename] |= bitarray('00000001')
                            if lig_atom.MatchesSMARTS('[Br]') and res_atom.IsHbondAcceptor() and dist <= 3.72:
                                for neighborHal in OBAtomAtomIter(lig_atom):
                                    if neighborHal.IsCarbon():
                                        angle = neighborHal.GetAngle(lig_atom, res_atom)
                                        if angle > 126.0:
                                            residue_dict[residuename] |= bitarray('00000001')
                            if lig_atom.MatchesSMARTS('[I]') and res_atom.IsHbondAcceptor() and dist <= 3.9:
                                for neighborHal in OBAtomAtomIter(lig_atom):
                                    if neighborHal.IsCarbon():
                                        angle = neighborHal.GetAngle(lig_atom, res_atom)
                                        if angle > 126.0:
                                            residue_dict[residuename] |= bitarray('00000001')
                        k += 1
                j += 1

    print 'done with other interactions with conformation: '+str(conform+1)

    output_holder.append([])
    for residue in resname:
        output_holder[conform].append(str(residue_dict[residue]).split("'")[1])


print 'done with interaction analysis'

output_arrays = []
output_res = []
remove_res = []
#Go through all of the fingerprint to see if any residues don't have any inteactions
for i in range(len(output_holder[0])):
    remove_res.append(True)
    for j in range(len(output_holder)):
        if '1' in output_holder[j][i]:
            remove_res[i] = False
            output_res.append(resname[i])
#Remove the residues without any interactions from the fingerprint
for i in range(len(output_holder)):
    output_arrays.append([])
    for j in range(len(remove_res)):
        if not remove_res[j]:
            output_arrays[i].append(output_holder[i][j])

output_res = list(OrderedDict.fromkeys(output_res))
graph_xlabels = []
for i in range(len(output_res)):
    graph_xlabels.append(output_res[i][3:])
graph_array = []
for i in range(len(output_arrays)):
    graph_array.append([])
    for j in range(len(output_arrays[0])):
        for k in range(len(output_arrays[0][0])):
            graph_array[i].append(int(output_arrays[i][j][k]))

#Write to file
output_file = open(output, 'w')
for i in range(len(output_res)):
    output_res[i] = output_res[i].ljust(9)
    output_file.write(output_res[i])
output_file.write('\n')
for i in range(len(output_arrays)):
    output_file.write('|'.join(output_arrays[i]))
    output_file.write('\n')
output_file.close()

#Create graph
cmap = mpl.colors.ListedColormap(['white','red'])
bounds=[0,0.5,1]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
max_axis = max(len(graph_array[0]), len(graph_array))
min_axis = min(len(graph_array[0]), len(graph_array))
ratio = int(max_axis / min_axis)
aspect_ratio = int(ratio / 3)
if aspect_ratio == 0:
    fig = pyplot.figure(figsize=(7*ratio,7))
    aspect_ratio = 'auto'
else:
    fig = pyplot.figure(figsize=(7*aspect_ratio,7))

ax = pyplot.axes(yticks=[])
majorLocator   = FixedLocator(np.linspace(0,len(graph_array[0]),len(output_res)+1))
ax.xaxis.set_major_locator(majorLocator)
img = pyplot.imshow(graph_array,interpolation='nearest',aspect=aspect_ratio,
                    cmap = cmap,norm=norm)
pyplot.tick_params(axis='x', which='both', bottom='off', top='off')
for i in range(len(output_arrays)):
    pyplot.axhline(y=i+0.5, color='k')
for i in range(len(graph_array[0])):
    if i%8 == 0:
        pyplot.axvline(x=i-0.5, color='k', linestyle='-', alpha=0.2)
ax.xaxis.set_ticklabels(graph_xlabels, rotation='vertical')
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_horizontalalignment('left')
if not graph_title == "":
    pyplot.title(graph_title)

pyplot.savefig(graph_name, format='tiff', dpi=300)
pyplot.show()







