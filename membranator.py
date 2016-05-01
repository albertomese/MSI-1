import os
import sys
import re
import argparse
import itertools
import logging
import time
from Bio.PDB import *
from Bio.PDB.MMCIF2Dict import *
from htmd import *
from numpy import *



def get_secondary_structure(dire, filename, mol):

    """
    Here we will select secondary structures suitable for being transmembrane domains in a list of lists
    """
    tm_helix = list()
    tm_sheet = list()
    tm_type = str()
    hidro_helix = 0
    hidro_sheet = 0
    hidropho_scale = {"ASP": -1.23, "GLU": -2.02, "ASN": -0.42, "GLN": -0.58, "LYS": -0.99, "ARG": -0.81, "HIS": -0.96, "GLY": -0.01, "PRO": -0.45, "SER": -0.13, "THR": -0.14, "CYS":  0.24, "MET":  0.23, "MSE":  0.23, "ALA": -0.17, "VAL": -0.07, "ILE":  0.31, "LEU":  0.56, "PHE":  1.13, "TRP":  1.85, "TYR":  0.94}
    sys.stderr.write("Computing secondary structure of " + filename + " \n")
    try:
        try:
            mmcif_dict = MMCIF2Dict(dire + filename + ".cif")
            tm_info = [mmcif_dict['_struct_conf.beg_auth_asym_id'], mmcif_dict['_struct_conf.beg_auth_seq_id'], mmcif_dict['_struct_conf.end_auth_seq_id']]
            ###Genero listas con la cadena, el inicio y el final de cada helice, luego las proceso para seleccionar las mas largas de 20
            ###tm_info es una lista que a su vez contiene listas -> tm_info[1] contiene cadena, inicio y fin para la segunda 
            i = 0
            while i < len(tm_info[0]):
                if int(tm_info[2][i]) - int(tm_info[1][i]) >= 14:
                    l = [tm_info[0][i], tm_info[1][i], tm_info[2][i]]
                    tm_helix.append(l)
                    ###Ahora tenemos una lista de tuples que contiene las coordenadas de nuestros dominios transmembrana -> funciona bien
                i += 1
        except:
            sys.stderr.write("The protein " + filename + " has no alpha helices. \n")


        try:
            tm_info = [mmcif_dict['_struct_sheet_range.beg_auth_asym_id'], mmcif_dict['_struct_sheet_range.beg_auth_seq_id'], mmcif_dict['_struct_sheet_range.end_auth_seq_id']]
            ###Genero listas con la cadena, el inicio y el final de cada helice, luego las proceso para seleccionar las mas largas de 20
            ###tm_info es una lista que a su vez contiene listas -> tm_info[1] contiene cadena, inicio y fin para la segunda 
            i = 0
            while i < len(tm_info[0]):
                if int(tm_info[2][i]) - int(tm_info[1][i]) >= 8:
                    l = [tm_info[0][i], tm_info[1][i], tm_info[2][i]]
                    tm_sheet.append(l)
                i += 1
        except:
            sys.stderr.write("The protein " + filename + " has no beta sheets. \n")

    except:
        sys.stderr.write("Secondary structure of " + filename + " failed. \n")

    ###Now we have identified the alpha helix and the beta sheet domains in our protein, here we will find if the transmembrane domains of our 
    ###protein are beta sheet or alpha helical.
    if tm_helix == []:
        tm_type = "sheet"

    elif tm_sheet == []:
        tm_type = "helix"

    elif tm_helix == [] and tm_sheet == []:
        sys.stderr.write("Secondary structure of " + filename + " failed. \n")
        return

    else:
        sys.stderr.write("Computing hidrophobicity of " + filename + " \n")
        i = 0
        for tm in tm_helix:
            for aa in range (int(tm[1]), int(tm[2]), 1):
                res = res = mol.get("resname", "name CA and chain " + str(tm[0]) + " and resid " + str(aa) + " and protein")
                val = hidropho_scale[res[0]]
                hidro_helix += float(val)
                i += 1
        hidro_helix = hidro_helix/i
    

        i = 0
        for tm in tm_sheet: ####Substituir las means por comandos de numpy
            for aa in range (int(tm[1]), int(tm[2]), 1):
                res = mol.get("resname", "name CA and chain " + str(tm[0]) + " and resid " + str(aa) + " and protein")
                val = hidropho_scale[res[0]]
                hidro_sheet += float(val)
                i += 1
        hidro_sheet = hidro_sheet/i
        

        if hidro_helix >= hidro_sheet:
            tm_type = "helix"
        else:
            tm_type = "sheet"

    if tm_type == "helix":
        sec_str = tm_helix
    else:
        sec_str = tm_sheet

    #here I process the list of secondary structures: removing repeated elements and ordering them 
    sec_str.sort()
    sec_str = list(sec_str for sec_str,_ in itertools.groupby(sec_str))
    sec_str.sort(key=lambda x: int(x[1]))

    tup = (tm_type, sec_str)
    return tup

def get_helix_vector(mol, sec_str):
    """
    Function to obtain the perpendicular vector to the membrane by using the C=O of the
    aminoacids in the helix
    """
    vectors=[]
    whole_vect=[]
    i = 1
    for n in range(len(sec_str)):
        for a in range(int(sec_str[n][1])+1, int(sec_str[n][2])-1):
            c_atoms = mol.get("coords","name C and chain " + str(sec_str[n][0]) + " and resid " + str(a) + " and protein")
            if len(c_atoms) != 3:
                try:
                    c_atoms = c_atoms[0]
                except:
                    continue 
            o_atoms = mol.get("coords","name O and chain " + str(sec_str[n][0]) + " and resid " + str(a) + " and protein")
            if len(o_atoms) != 3:
                try:
                    o_atoms = o_atoms[0]
                except:
                    continue 

            vectors.append(o_atoms - c_atoms)
        helix_vector=np.mean(vectors, axis=0)
        #As contiguous helices may point in opposite directions, we invert the vectors corresponding with a pair transmembrane domain 
        if i % 2 == 0:
            helix_vector = helix_vector * [-1, -1, -1]
        i += 1
        whole_vect.append(helix_vector)
    final_vect = np.array(np.mean(whole_vect, axis=0))
    return final_vect


def get_sheet_vector_CA_CA_version(mol, sec_str):
    """
    Here we will compute the perpendicular vector to the membrane for beta sheet proteins by computing the average vector between distal
    residues in each beta sheet.
    """
    #We start by selecting CA located at the end of each beta sheet segment, neglecting the two last residues from each sheet 
    vector_list = list()
    i = 1
    for element in sec_str:
        start = mol.get("coords", "name CA and resid " + str(int(element[1])+2) + " and protein")
        if len(start) != 3:
                try:
                    start = start[0]
                except:
                    continue 
        end = mol.get("coords", "name CA and resid " + str(int(element[2])-2) + " and protein")
        if len(end) != 3:
                try:
                    end = end[0]
                except:
                    continue 
        vec = end - start
        #As beta barrels have contiguous sheets in oposite directions, we will compute the oposite vector for even sheets, so all vectors point in the same direction
        if i % 2 == 0:
            vec = vec * [-1, -1, -1]
        #Now, all vectors must be unitary for obtaining an unbiased mean
        norm = (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
        vec[0] = vec[0]/norm
        vec[1] = vec[1]/norm
        vec[2] = vec[2]/norm
        vector_list.append(vec)
        i += 1

    sheet_vector = np.mean(vector_list, axis=0)
    return sheet_vector
 

def get_sheet_vector_consecutive_CA_version(mol, sec_str):
    """
    Here we will compute the perpendicular vector to the membrane for beta sheet proteins by computing the mean vector between CA 
    of adjacent residues in the transmembrane domains.
    """
    #We start by selecting CA located at the end of each beta sheet segment, neglecting the last residues from each sheet 
    vector_list = list()
    super_vector_list = list()
    i = 1
    for element in sec_str:
        for aa in range((int(element[1])+3), (int(element[2])-3), 1):
            c1 = mol.get("coords", "name CA and chain " + element[0] + " and resid " + str(aa) + " and protein")
            if len(c1) != 3:
                try:
                    c1 = c1[0]
                except:
                    continue 
            c2 = mol.get("coords", "name CA and chain " + element[0] + " and resid " + str(aa+1) + " and protein")
            if len(c2) != 3:
                try:
                    c2 = c2[0]
                except:
                    continue 
            vec = c2 - c1
            vector_list.append(vec)

        super_vector = np.mean(vector_list, axis=0)
        norm = (super_vector[0]**2 + super_vector[1]**2 + super_vector[2]**2)**0.5
        super_vector[0] = super_vector[0]/norm
        super_vector[1] = super_vector[1]/norm
        super_vector[2] = super_vector[2]/norm
        #As beta barrels have contiguous sheets in oposite directions, we will compute the oposite vector for even sheets, so all vectors point in the same direction
        if i % 2 == 0:
            super_vector = super_vector * [-1, -1, -1]
        i += 1
        super_vector_list.append(super_vector)
        vector_list = list()
    final_vector = np.mean(super_vector_list, axis=0)
    return final_vector

def get_sheet_vector_inertia_tensor_version(mol, sec_str):
    """
    Here we will compute the perpendicular vector to the membrane for beta sheet proteins by computing the inertia tensor matrix.
    For that we asume the beta barrel as a cylinder, and we consider the CA of the beta sheets for determining mass distribution 
    in this cylinder.
    """
    #We start by selecting CA located at the end of each beta sheet segment, neglecting the two last residues from each sheet 
    coordinates = list()
    I = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for element in sec_str:
        for aa in range(int(element[1]), int(element[2]), 1):
            try:
                res_coord = mol.get("coords", "name CA and chain " + element[0] + " and resid " + str(aa) + " and protein")
                if len(res_coord) != 3:
                    try:
                        res_coord = res_coord[0]
                    except:
                        continue 
                coordinates.append(res_coord)
            except:
                continue
    
    for element in coordinates: 

        #We compute the matrix for the inertia tensor
        I[0][0] += element[1]**2 + element[2]**2
        I[0][1] -= element[0]*element[1]
        I[0][2] -= element[0]*element[2]
        I[1][0] -= element[1]*element[0]
        I[1][1] += element[0]**2 + element[2]**2
        I[1][2] -= element[1]*element[2]
        I[2][0] -= element[2]*element[0]
        I[2][1] -= element[2]*element[1]
        I[2][2] += element[1]**2 + element[0]**2

    diagonalization = np.linalg.eig(I)
    
    eigenvals = list(diagonalization[0])
    eigenvalue = min(diagonalization[0])
    eigencoord = eigenvals.index(eigenvalue)
    #The output of the diagonalization function has two main arrays, the first containing the eigenvalues and the second 
    #containing the eigenvectors. So we select the second array and on it we get the eigenvectors corresponding to the lowest eigenvalue
    final_vector = diagonalization[1][eigencoord]
    return final_vector

def get_sheet_vector_ring_version(mol, sec_str):
    """
    Here we will compute the perpendicular vector to the membrane for beta sheet proteins by computing the average perpendicular vector to the 
    aromatic rings in the transmembrane domains
    """
    vec = [0, 0, 0]
    vector_list = []
    vec_x = list()
    vec_y = list()
    vec_z = list()
    for element in sec_str:
        for aa in range(int(element[1]), int(element[2])-1, 1):
            if mol.get('resname', 'name CA and resid ' + str(aa) + ' and protein') == 'HIS':
                try:
                    cg = mol.get('coords', 'name CG and resid ' + str(aa) + ' and protein')
                    if len(cg) != 3:
                        try:
                            cg = cg[0]
                        except:
                            continue 
                    cd1 = mol.get('coords', 'name ND1 and resid ' + str(aa) + ' and protein')
                    if len(cd1) != 3:
                        try:
                            cd1 = cd1[0]
                        except:
                            continue
                    cd2 = mol.get('coords', 'name CD2 and resid ' + str(aa) + ' and protein')
                    if len(cd2) != 3:
                        try:
                            cd2 = cd2[0]
                        except:
                            continue
                    v1 = cd1 - cg
                    v2 = cd2 - cg
                    #Now we compute the vector product for v1 and v2
                    vec_x.append(v1[1]*v2[2] + v2[1]*v1[2])
                    vec_y.append(-(v1[0]*v2[2] + v2[0]*v1[2]))
                    vec_z.append(v1[0]*v2[1] + v2[0]*v1[1])
                except:
                    continue

            elif (mol.get('resname', 'name CA and resid ' + str(aa) + ' and protein') == 'TYR') or (mol.get('resname', 'name CA and resid ' + str(aa) + ' and protein') == 'TRP') or (mol.get('resname', 'name CA and resid ' + str(aa) + ' and protein') == 'PHE'):
                try:
                    cg = mol.get('coords', 'name CG and resid ' + str(aa) + ' and protein')
                    if len(cg) != 3:
                        try:
                            cg = cg[0]
                        except:
                            continue
                    cd1 = mol.get('coords', 'name CD1 and resid ' + str(aa) + ' and protein')
                    if len(cd1) != 3:
                        try:
                            cd1 = cd1[0]
                        except:
                            continue
                    cd2 = mol.get('coords', 'name CD2 and resid ' + str(aa) + ' and protein')
                    if len(cd2) != 3:
                        try:
                            cd2 = cd2[0]
                        except:
                            continue
                    v1 = cd1 - cg
                    v2 = cd2 - cg
                    #Now we compute the vector product for v1 and v2
                    vec_x.append(v1[1]*v2[2] + v2[1]*v1[2])
                    vec_y.append(-(v1[0]*v2[2] + v2[0]*v1[2]))
                    vec_z.append(v1[0]*v2[1] + v2[0]*v1[1])
                except:
                    continue

    x = sum(vec_x)/len(vec_x)
    y = sum(vec_y)/len(vec_y)
    z = sum(vec_z)/len(vec_z)
    final_vector = [x, y, z]

    
    return final_vector

def rotate_prot(vector, mol):
    """
    Function to rotate the protein using the 
    """
    rotate_np = np.cross(vector, np.array([[0], [0], [1]]), axis=0)
    radius = np.sqrt((vector[0]**2)+(vector[1]**2)+(vector[2]**2))
    theta = np.arccos(vector[2]/radius)
    rotate_vec = rotate_np.ravel()
    mol.rotate(rotate_vec,theta)

def hydrophobicity_calculator(minimum, maximum, mol, jump):
    aa_hidro=["ALA","LEU","VAL","ILE","TRP","PHE","PRO","MET"]
    counter = minimum
    counter2 = minimum + jump 
    hidro_list=[]
    hidro_counter = 0
    while counter <= maximum:
        res_list = mol.get('resname', 'z> ' + str(counter) + ' and z< ' + str(counter2) + ' and protein')
        for element in res_list:
            if element in aa_hidro:
                hidro_counter += 1
        hidro_list.append(hidro_counter)
        hidro_counter = 0
        counter += jump
        counter2 += jump
    density=[]
    x = 0
    for a in hidro_list:
        if x == 0:
            density.append((a+hidro_list[x+1])/sum(hidro_list))
        if x+1 == len(hidro_list):
            density.append((a+hidro_list[x-1])/sum(hidro_list))
        else:
            density.append((a+hidro_list[x+1]+hidro_list[x-1])/sum(hidro_list))
        x +=1
    mem = False
    for a in range(len(density)):
        if density[a] >= np.mean(density,axis=0) and not mem:
            mem = True 
            bot = a
        elif density[a] >= np.mean(density,axis=0) and mem == True:
            top = a
            if abs(top - bot) > 20:
                break
    
    bot = bot*2 + minimum
    top = top*2 + minimum

    return (hidro_list,top,bot)

def dummy_leaflet(pmin, pmax,z,spacing, mol):
    membrane=Molecule()
    Ni = int(((pmax[0] - pmin[0])/spacing)+1)
    Nj = int(((pmax[1] - pmin[1])/spacing)+1)
    membrane.empty(Ni*Nj)
    membrane.set('record', 'HETATM')
    membrane.set('name', 'MEM')
    membrane.set('resname', 'MRES')
    resids=np.array(range(Ni*Nj))
    membrane.set('resid', resids)
    coordinate_list = []
    for col in range(Ni):
        for row in range(Nj):
            coordinate_list.append([pmin[0]+(col*spacing),pmin[1]+(row*spacing),z])
    membrane.set('coords', coordinate_list)
    mol.append(membrane)
    

    return (membrane)


def membrane_builder(dire, filename, bmode):

    """
    Here we compute the plane corresponding to the center of the lipidic membrane
    """
    mol = Molecule(dire + filename + '.pdb') 
    structure = get_secondary_structure(dire, filename, mol)
    tm_type = structure[0]
    sec_str = structure[1]
    if tm_type == "helix":
        sys.stderr.write(filename + " has an alpha helical transmembrane domain \n")
        vector = get_helix_vector(mol, sec_str)
    else:
        sys.stderr.write(filename + " has an beta sheet transmembrane domain \n")
        if bmode == "0":
            sys.stderr.write("Calculating vector with the CA - CA method \n")
            vector = get_sheet_vector_CA_CA_version(mol, sec_str)
        elif bmode == "1":
            sys.stderr.write("Calculating vector with the consecutive CA method \n")
            vector = get_sheet_vector_consecutive_CA_version(mol, sec_str)
        elif bmode == "2":
            sys.stderr.write("Calculating vector with inertia tensor method \n")
            vector = get_sheet_vector_inertia_tensor_version(mol, sec_str)
        elif bmode == "3":
            sys.stderr.write("Calculating vector with the aromatic rings method \n")
            vector = get_sheet_vector_ring_version(mol, sec_str)

    #Here we load the rotation and the visualization
    sys.stderr.write("Rotating the protein... \n")
    rotate_prot(vector,mol)
    
    coords=mol.get("coords")
    xyz=list(zip(*coords))
    max_z = int(round(max(xyz[2])))
    min_z = int(min(xyz[2]))
    sys.stderr.write("Calculating hydrophobic areas\n")
    hidro_list,top,bot = hydrophobicity_calculator(min_z, max_z, mol, 2)
    sys.stderr.write("Generating the membranes... \n")
    max_x = int(round(max(xyz[0])))+5
    max_y = int(round(max(xyz[1])))+5
    min_x = int(min(xyz[0]))-5
    min_y = int(min(xyz[1]))-5
    dummy_leaflet([min_x,min_y],[max_x,max_y], bot, 2, mol)
    dummy_leaflet([min_x,min_y],[max_x,max_y], top, 2, mol)
    #mol.write('./transmembrane.pdb')
    #mol.view()
    #import time
    #time.sleep(60000)

    return mol


def argparse_creator():

    parser = argparse.ArgumentParser(description="""This program gets a pdb and a cif structure and returns a pdb file with the protein placed 
                                                    within two layers of dummy atoms corresponding to a lipidic membrane.""")

    parser.add_argument('-i','--input',
                        dest = "infile",
                        action = "store",
                        default = None,
                        help = "Input file or directory. Files must have only one chain")

    parser.add_argument('-o','--output',
                        dest = "outfile",
                        action = "store",
                        default = "output",
                        help = "Path and output files.")

    parser.add_argument('-b','--beta',
                        dest = "bmode",
                        action = "store",
                        default = "0",
                        help = "Select the algorithm you prefer for the orientation of the protein, if it has a beta sheet transmembrane domain.")

    options = parser.parse_args()

    infile = options.infile
    output = options.outfile
    bmode = options.bmode

    arg = (infile, output, bmode)
    return arg

if __name__ == '__main__': ###esto permite que pueda ejecutar o extraer funciones de este modulo sin necesidad de ejecutar todas las funciones del modulo

    argv = argparse_creator()
    if len(argv) >= 2:
        dire = argv[0] ###esta variable contendrá el path
        out = argv[1]
        bmode = argv[2]
        filenames = list()
        if os.path.isdir(dire): ###devuelve true si el contenido es un directorio
            b = re.compile(".*\.cif")
            files = list(os.listdir(str(dire)))
            for element in files:
                if b.match(str(element)):
                    filenames.append(element[:-4])
            files_found = len(filenames)
            sys.stderr.write(str(files_found) + " files found. \n")
            i = 1
            for filename in filenames:
                sys.stderr.write("Working with " + filename + ". \n")
                mol = membrane_builder(dire, filename, bmode) 
                mol.write('./' + str(out) + str(i) + '.pdb')
                i += 1
                

        else:
            raise ValueError("You are missing the directory containing the database")
    else:
        raise ValueError("You are missing the directory containing the database")