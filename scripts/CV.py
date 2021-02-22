#!/usr/bin/env python3

import math


# CV calculation
def distance(atom1, atom2):
    x = atom1[0] - atom2[0]
    y = atom1[1] - atom2[1]
    z = atom1[2] - atom2[2]

    return math.sqrt(x ** 2 + y ** 2 + z ** 2)


# calcule la liste d'atomes dans un rayon de rc de l'atome de reference de coords (x,y,z)
def Env_i(xref, yref, zref, rc, atomlist):
    """input:
            - xi, yi, zi  --> coords of atomi (floats)
            - rc --> radius in angstm (float)
            - atomilst --> list of coords of all atoms (list)
        output:
            - nei --> list of coords of all atomi's neighbors (list)
    """
    nei = []
    for i in range(len(atomlist)):
        if distance((xref, yref, zref), atomlist[i]) < rc:
            nei.append(atomlist[i])
    return nei


# calcule le vecteur defini entre 2 atomes
def Rij(atom1, atom2):
    return (atom1[0] - atom2[0], atom1[1] - atom2[1], atom1[2] - atom2[2])


# calcule la norme du vecteur de coords vx, vy,vz
def Norme(vx, vy, vz):
    return math.sqrt(vx ** 2 + vy ** 2 + vz ** 2)


# calcule CV de l'atome i
def CVi(xi, yi, zi, rc, atomlist):
    """input:
            - xi, yi, zi  --> coords of atomi (floats)
            - rc --> radius in angstm (float)
            - atomilst --> list of coords of all atoms (list)
        output:
            - CV --> CV of atomi (float)
    """

    neighbors = Env_i(xi, yi, zi, rc, atomlist)
    SVX = 0.0
    SVY = 0.0
    SVZ = 0.0

    for j in range(len(neighbors)):
        vx, vy, vz = Rij((xi, yi, zi), neighbors[j])
        norm = Norme(vx, vy, vz)
        if norm != 0:
            SVX += vx / norm
            SVY += vy / norm
            SVZ += vz / norm

    CV = 1 - Norme(SVX, SVY, SVZ) / len(neighbors)
    if Norme(SVX, SVY, SVZ) / len(neighbors) > 1:
        print("CV < 0 !")
        print("Norme: ", Norme(SVX, SVY, SVZ))
        print("nb atomes voisins : ", len(neighbors), "\n")
    return CV


# calcule la CV de chaque atome
def CV_AllRes(atomlist, dPDB, chainID, rc):
    """
        input:
            - atomilst --> list of coords of all atoms (list)
            - dPDB --> dico PDB (dict)
            - chainID --> name of the chain to treat (string)
            - rc --> radius in angstm (float)
         output:
             - returns nothing. All CVs are stored in the "bfactor" key of dPDB
    """

    for resi in dPDB[chainID]["reslist"]:
        CV_resi = 0

        for atomi in dPDB[chainID][resi]["atomlist"]:
            # recupere coords de l'atome i (pas utile, on aurait pu donner directement a la fct CVi dPDB[chainID][resi][atomi]["x"] etc mais code plus lisible)
            xi = dPDB[chainID][resi][atomi]["x"]
            yi = dPDB[chainID][resi][atomi]["y"]
            zi = dPDB[chainID][resi][atomi]["z"]

            # calcule la CV de chq atomi du resi
            CV_atomi = CVi(xi, yi, zi, rc, atomlist)
            CV_resi += CV_atomi

        # calcule la moyenne des CV des atomes du resi
        CV_resi = CV_resi / len(dPDB[chainID][resi]["atomlist"])

        for atomi in dPDB[chainID][resi]["atomlist"]:
            dPDB[chainID][resi][atomi]["bfactor"] = CV_resi  # utile pour ecrire ensuite la CV moyenne dans le PDB (on la transfere a chq atome du resi en question

