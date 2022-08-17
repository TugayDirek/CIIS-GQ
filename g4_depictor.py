import glob
import os
import time
import numpy as np
from copy import deepcopy
import turtle
# import image
from PIL import Image
import PIL.Image
import svgwrite
# from svg_turtle import SvgTurtle
# import pyvips
from PIL import Image
# import pillow

pwdict = ["GLU", "AGLU","BGLU","CGLU","DGLU",
"LEU", "ALEU","BLEU","CLEU","DLEU",
"THR", "ATHR","BTHR","CTHR","DTHR",
"PRO", "APRO","BPRO","CPRO","DPRO",
"ASP", "AASP","BASP","CASP","DASP",
"GLN", "AGLN","BGLN","CGLN","DGLN",
"HIS", "AHIS","BHIS","CHIS","DHIS",
"PHE", "APHE","BPHE","CPHE","DPHE",
"ILE", "AILE","BILE","CILE","DILE",
"MET", "AMET","BMET","CMET","DMET",
"SER", "ASER","BSER","CSER","DSER",
"TYR", "ATYR","BTYR","CTYR","DTYR",
"ASN", "AASN","BASN","CASN","DASN",
"LYS", "ALYS","BLYS","CLYS","DLYS",
"ARG", "AARG","BARG","CARG","DARG",
"ALA", "AALA","BALA","CALA","DALA",
"VAL", "AVAL","BVAL","CVAL","DVAL",
"GLY", "AGLY","BGLY","CGLY","DGLY",
"TRP", "ATRP","BTRP","CTRP","DTRP",
"CYS", "ACYS","BCYS","CCYS","DCYS",
]

modifiefGList = ["BGM", "8OG", "GFL", "GF2", "GF0", "LCG", "0G" ]

def drawStructure(tetrad_and_loop, pdb_name,ONZ,tetrad_onz, atoms_sequence):
    key = tetrad_and_loop[pdb_name]
    sequence = atoms_sequence[pdb_name]
    tetrads = key[0]
    loop = key[1]
    global_G=10000
    bulge_bases = []




    for index in range(0, len(tetrads)):
        for i in range(0, len(tetrads[index])):
            tetrads[index][i] = int(tetrads[index][i])

    foundG4s = {}
    first_index = -1

    # assign each guanin in tetrad to the index it resides in to a map
    for g4 in tetrads:
        first_index += 1

        for point in range(0, len(g4)):
            foundG4s[g4[point]] = first_index

    pen = turtle.Turtle()


    screen = turtle.Screen()

    pen.hideturtle()
    screen.title(pdb_name)
    style = ("Courier", 13, "bold")
    style2 = ("Courier", 11, "italic")
    pen.penup()
    pen.color('brown')
    counterONZ = 0
    for tetradType in ONZ[pdb_name]:
        pen.setpos(-250, 50+counterONZ*60)

        pen.write(tetradType+")", font=style2, align="right")
        counterONZ+=1


    pen.setpos(-50, 0)
    tetrad_list = ""
    for tetrads in tetrad_onz[pdb_name]:
        tetrad_list = tetrad_list + "["
        for tetrad in tetrads:
            tetrad_list = tetrad_list + str(tetrad) + ","
        tetrad_list = tetrad_list[:-1]
        tetrad_list = tetrad_list + "],"

    tetrad_list = tetrad_list[:-1]


    pen.write("Tetrads: "+tetrad_list, font=style2, align="center")
    pen.setpos(-200, 50)
    #pen.pendown()
    pen.color('black')
    y_positions = []
    y_positions.append(50)
    y_positions.append(110)
    y_positions.append(170)
    y_positions.append(230)
    y_positions.append(300)


    enter = False
    loop_type = ""
    last_x = 0

    for element in loop:


        x = -200
        for element2 in element:

            if (element2.isnumeric()):

                g_number = int(element2)


                if (abs(global_G - g_number) > 1 and abs(global_G - g_number) < 800 and not enter):
                    if(foundG4s[int(global_G)] - foundG4s[int(g_number)]<0):
                        pen.setpos(x, ((y_positions[foundG4s[int(global_G)]] + y_positions[foundG4s[int(g_number)]]) / 2) - 7)
                        pen.penup()

                        pen2 = turtle.Turtle()
                        pen2.hideturtle()
                        pen2.penup()
                        pen2.setpos(x, ((y_positions[foundG4s[int(global_G)]] + y_positions[foundG4s[int(g_number)]]) / 2) - 7)
                        pen2.pendown()
                        pen2.circle(7, 180)

                        pen.setpos(x, ((y_positions[foundG4s[int(global_G)]] + y_positions[foundG4s[int(g_number)]]) / 2) + 7)
                        pen.pendown()
                        base_numbers = []
                        for b in range(global_G + 1, g_number):
                            base_numbers.append(b)
                        bulge_bases.append(base_numbers)
                    else:
                        pen.setpos(x, ((y_positions[foundG4s[int(global_G)]] + y_positions[foundG4s[int(g_number)]]) / 2) + 7)
                        pen.penup()

                        pen2 = turtle.Turtle()
                        pen2.hideturtle()
                        pen2.penup()
                        pen2.setpos(x, ((y_positions[foundG4s[int(global_G)]] + y_positions[foundG4s[int(g_number)]]) / 2) + 7)
                        pen2.pendown()
                        pen2.circle(-7, 180)

                        pen.setpos(x, ((y_positions[foundG4s[int(global_G)]] + y_positions[foundG4s[int(g_number)]]) / 2) - 7)
                        pen.pendown()
                        base_numbers = []
                        for b in range(global_G+1, g_number):
                            base_numbers.append(b)
                        bulge_bases.append(base_numbers)



                if (enter):
                    z = last_x + 55


                    if (loop_type.__contains__("reversal")):
                        pen.color("red")
                    else:
                        if (loop_type.__contains__("lateral")):
                            pen.color("green")
                        else:
                            pen.color("blue")


                    pen.setpos(z, (y_positions[foundG4s[int(element2)]] + y_positions[foundG4s[int(global_G)]]) / 2)
                    pen.write(loop_type, align="center")
                    enter = False
                    pen.color("black")

                pen.setpos(x, y_positions[foundG4s[int(element2)]])



            else:
                last_x = x
                x += 110 # go towards right after each loop
                enter = True
                loop_type = element2
                continue
            pen.pendown()
            pen.write(g_number, font=("Courier", 10, "bold"))
            letter = sequence[str(g_number)]
            if(sequence[str(g_number)].__contains__("GM")):
                letter="BGM"
            elif (sequence[str(g_number)].__contains__("FL")):
                letter = "GFL"
            elif (sequence[str(g_number)].__contains__("F0")):
                letter = "GF0"
            elif (sequence[str(g_number)].__contains__("OG")):
                letter = "8OG"
            elif (sequence[str(g_number)].__contains__("F2")):
                letter = "GF2"
            elif (sequence[str(g_number)].__contains__("CG")):
                letter = "LCG"

            pen.write(letter, font=style, align="right")


            global_G=g_number


    pen.penup()
    pen.setpos(-40, -30)

    bulges = ""
    for bases in bulge_bases:
        bulges = bulges + "["
        for base in bases:
            if(str(base) in sequence):
                bulges = bulges + str(base) + sequence[str(base)] + ","
            else:
                bulges = bulges + str(base) + ","

        bulges = bulges[:-1]
        bulges = bulges + "],"

    bulges = bulges[:-1]

    if(len(bulge_bases)>0):
        pen.color("brown")
        pen.write("Bulge bases: "+ bulges, font=style2, align="center")
        print("Bulge Bases: ",bulge_bases)


    turtle.delay(100)
    #time.sleep(1)
    # pen.onclick(time.sleep(5))
    #pen.clear()
    turtle.done()


# return 3D coordinates of guanines and sequence of the DNA/RNA
def readFromFiles(files):
    counter = 0
    counterRNA = 0
    counterDNA = 0
    typeOfNucleicAcid = ""
    atomsDNA = {}
    atomsRNA = {}
    atomsDNA_sequence = {}
    atomsRNA_sequence = {}

    for filename in files:

        pdbName = filename.split("\\")[-1].split(".")[0]#[1]
        with open(filename, 'r') as f:

            cont = True
            first = ""
            coordinatesDict = {}
            check = ""
            counter += 1
            ligand_counter = 0

            for line in f:
                '''''
                '''
                if (line.split()[0].__contains__("HEADER") and line.__contains__("DNA")):

                    atomsDNA[pdbName] = {}
                    typeOfNucleicAcid = "DNA"
                    counterDNA += 1

                    atomsDNA_sequence[pdbName] = {}


                if (line.split()[0].__contains__("HEADER") and line.__contains__("RNA")):
                    atomsRNA[pdbName] = {}
                    typeOfNucleicAcid = "RNA"

                    counterRNA += 1

                    atomsRNA_sequence[pdbName] = {}

                if (line.split()[0].__contains__("HEADER") and (line.__contains__("DNA") and line.__contains__("RNA"))):
                    typeOfNucleicAcid = "DNA/RNA"

                coordinates = []
                # read until the first atom is seen
                if ((line.split()[0].__contains__("ATOM") or line.split()[0].__contains__("HETATM")) and line.split()[
                    1].__contains__("1")):
                    cont = False

                if (cont):
                    continue

                # finish reading after seeing the end of the protein
                if line.__contains__("ENDMDL") or line.__contains__("END"):
                    break

                # print(filename)
                list = line.split()

                # if there is another structure we don't want, delete it
                if ((check == ("TER") and (list[0] == "ATOM" and list[3] not in pwdict) )):  #    and len(list[3])!=3 or list[0] == "HETATM"  --  or (check == "ATOM" and list[0] == "HETATM")):  # (list[0].__contains__("ATOM") and str(line[17:20]) in pwdict) or
                    if (any((pdbName) in x for x in atomsRNA.keys())):
                        del (atomsRNA[pdbName])
                    if (any((pdbName) in x for x in atomsDNA.keys())):
                        del (atomsDNA[pdbName])
                    if (any((pdbName) in x for x in atomsRNA_sequence.keys())):
                        del (atomsRNA_sequence[pdbName])
                    if (any((pdbName) in x for x in atomsDNA_sequence.keys())):
                        del (atomsDNA_sequence[pdbName])
                    break


                check = list[0]  # looks fot TER keyword, if it is present it means it is not monomer


                if (list[0].__contains__("HETATM")):
                    ligand_counter += 1


                if (ligand_counter > 5 and (
                        (pdbName in atomsDNA and typeOfNucleicAcid.__contains__("DNA")) or (
                        pdbName in atomsRNA and typeOfNucleicAcid.__contains__("RNA")))):

                    if (typeOfNucleicAcid.__eq__("DNA")):
                        atomsDNA_sequence.get(pdbName)[line[23:27].strip() + "_ligand"] = "ligand"
                    elif (typeOfNucleicAcid.__eq__("RNA")):
                        atomsRNA_sequence.get(pdbName)[line[23:27].strip() + "_ligand"] = "ligand"
                    elif (typeOfNucleicAcid.__contains__("DNA") and typeOfNucleicAcid.__contains__("RNA")):
                        atomsDNA_sequence.get(pdbName)[line[23:27].strip() + "_ligand"] = "ligand"
                        atomsRNA_sequence.get(pdbName)[line[23:27].strip() + "_ligand"] = "ligand"

                    ligand_counter = -1000

                # checks if the Guanin atoms (N1,N2,N7,O6) we looked inserted to coordinatedDict or not, if it is the case, continue until a new atom is seen
                if ((line[18:20].__contains__("G") or line[18:20].__contains__("FL") or line[18:20].strip().__eq__(
                        "F2") or line[17:20].strip().__eq__("GF0"))
                        and any(("N1") in x for x in coordinatesDict.keys()) and any(
                            ("N2") in x for x in coordinatesDict.keys()) and
                        any(("N7") in x for x in coordinatesDict.keys()) and any(
                            ("O6") in x for x in coordinatesDict.keys()) and first == line[23:27].strip()):
                    continue
                elif (first != line[
                               23:27].strip()):  # if the previous and the current atoms number is not equal, empty the distionary
                    coordinatesDict = {}
                    ligand_counter = 0

                first = line[23:27].strip()

                enter_coord = True  # to prevent taking coordinates two times in commons bases

                # HERE insert the coordinates of the atoms we searched for
                if ((list[0].__contains__("ATOM") or list[0].__contains__("HETATM")) and (
                        line[18:20].strip().__eq__("DG") or line[18:20].strip().__eq__("G") or line[
                                                                                               18:20].strip().__eq__(
                        "DT") or line[18:20].strip().__eq__(
                        "DA") or line[18:20].strip().__eq__("DC") or line[18:20].strip().__eq__("GM") or line[
                                                                                                         18:20].strip().__eq__(
                    "FL") or line[18:20].strip().__eq__("OG") or line[18:20].strip().__eq__("F2")
                        or line[17:20].strip().__eq__("GF0") or line[18:20].strip().__eq__(
                    "CG") or line[18:20].strip().__eq__(
                    "0G")) and typeOfNucleicAcid.__contains__("DNA") and atomsDNA.__contains__(pdbName)):

                    atomsDNA_sequence.get(pdbName)[line[23:27].strip()] = line[18:20].strip()

                    atomsDNA.get(pdbName)[line[23:27].strip() + "," + line[18:20].strip()] = line[18:20].strip()
                    if ((line[18:20].strip().__eq__("DG") or line[18:20].strip().__eq__("G") or line[
                                                                                                18:20].strip().__eq__(
                            "GM") or line[18:20].strip().__eq__("FL") or line[18:20].strip().__eq__("OG")
                         or line[18:20].strip().__eq__("F2") or line[17:20].strip().__eq__("GF0") or line[
                                                                                                     18:20].strip().__eq__(
                                "CG") or line[18:20].strip().__eq__("0G")) and (
                            line[13:15].strip() == ("N1") or line[13:15].strip() == ("N2") or line[13:15].strip() == (
                    "N7") or line[13:15].strip() == ("O6"))):
                        coordinates.append((line[31:38].strip()))
                        coordinates.append((line[39:46].strip()))
                        coordinates.append((line[47:56].strip()))
                        coordinatesDict[line[13:15].strip()] = coordinates
                        atomsDNA.get(pdbName)[line[23:27].strip() + "," + line[18:20].strip()] = coordinatesDict
                        enter_coord = False

                if ((list[0].__contains__("ATOM") or list[0].__contains__("HETATM")) and (
                        line[18:20].strip().__eq__("G") or line[18:20].strip().__eq__("DG") or line[
                                                                                               18:20].strip().__eq__(
                        "U") or
                        line[18:20].strip().__eq__("A") or line[18:20].strip().__eq__("C") or line[
                                                                                              18:20].strip().__eq__(
                    "GM")
                        or line[18:20].strip().__eq__("FL") or line[18:20].strip().__eq__("OG") or line[
                                                                                                   18:20].strip().__eq__(
                    "F2") or line[17:20].strip().__eq__("GF0") or line[18:20].strip().__eq__("CG") or line[18:20].strip().__eq__(
                    "0G")
                ) and typeOfNucleicAcid.__contains__("RNA") and atomsRNA.__contains__(pdbName)):

                    atomsRNA_sequence.get(pdbName)[line[23:27].strip()] = line[18:20].strip()
                    atomsRNA.get(pdbName)[line[23:27].strip() + "," + line[18:20].strip()] = line[18:20].strip()

                    if ((line[18:20].strip().strip().__eq__("G") or line[18:20].strip().__eq__("DG") or line[
                                                                                                        18:20].strip().__eq__(
                            "GM") or
                         line[18:20].strip().__eq__("FL") or line[18:20].strip().__eq__("OG") or line[
                                                                                                 18:20].strip().__eq__(
                                "F2") or line[17:20].strip().__eq__("GF0") or line[18:20].strip().__eq__("CG") or line[18:20].strip().__eq__(
                    "0G")) and (
                            line[13:15].strip() == ("N1") or line[13:15].strip() == ("N2") or line[13:15].strip() == (
                            "N7") or
                            line[13:15].strip() == ("O6"))):
                        # print(pdbName)
                        if (enter_coord):
                            coordinates.append(float(line[31:38].strip()))
                            coordinates.append(float(line[39:46].strip()))
                            coordinates.append(float(line[47:56].strip()))
                            coordinatesDict[line[13:15].strip()] = coordinates

                        atomsRNA.get(pdbName)[line[23:27].strip() + "," + line[18:20].strip()] = coordinatesDict

    # print("There are " + str(counterDNA) + " DNA, " + str(counterRNA) + " RNA files.")
    return atomsDNA, atomsRNA, atomsDNA_sequence, atomsRNA_sequence


# delete element with no data
def deletePoorData(atoms):
    for key in list(atoms.keys()):
        if (len(atoms.get(key)) == 0):
            del (atoms[key])


def calculateDistances_AngleMethod(atoms, max_bond_distance):
    distanceDict = {}
    for keys in atoms.keys():
        distanceDict[keys] = {}

        for key in atoms[keys]:
            if key.__contains__("G"):

                a = np.array([float(atoms[keys][key]["N1"][0]), float(atoms[keys][key]["N1"][1]),
                              float(atoms[keys][key]["N1"][2])])
                b = np.array([float(atoms[keys][key]["N2"][0]), float(atoms[keys][key]["N2"][1]),
                              float(atoms[keys][key]["N2"][2])])
                c = np.array([float(atoms[keys][key]["N7"][0]), float(atoms[keys][key]["N7"][1]),
                              float(atoms[keys][key]["N7"][2])])
                d = np.array([float(atoms[keys][key]["O6"][0]), float(atoms[keys][key]["O6"][1]),
                              float(atoms[keys][key]["O6"][2])])



                for key2 in atoms[keys]:
                    if key2.__contains__("G"):

                        a2 = np.array([float(atoms[keys][key2]["N1"][0]), float(atoms[keys][key2]["N1"][1]),
                                       float(atoms[keys][key2]["N1"][2])])
                        b2 = np.array([float(atoms[keys][key2]["N2"][0]), float(atoms[keys][key2]["N2"][1]),
                                       float(atoms[keys][key2]["N2"][2])])
                        c2 = np.array([float(atoms[keys][key2]["N7"][0]), float(atoms[keys][key2]["N7"][1]),
                                       float(atoms[keys][key2]["N7"][2])])
                        d2 = np.array([float(atoms[keys][key2]["O6"][0]), float(atoms[keys][key2]["O6"][1]),
                                       float(atoms[keys][key2]["O6"][2])])

                        if (key is not key2):

                            N2_N7 = b - c2
                            N1_O6 = a - d2

                            N7_N2 = c - b2
                            O6_N1 = d - a2




                            if (np.linalg.norm(N2_N7) <=max_bond_distance and np.linalg.norm(N1_O6) <= max_bond_distance):
                                # print(np.linalg.norm(N2_N7))
                                # print(np.linalg.norm(N1_O6))
                                distanceDict.get(keys)["," + key + "," + key2 + ",N2_N7"] = N2_N7
                                distanceDict.get(keys)["," + key + "," + key2 + ",N1_O6"] = N1_O6

                            if (np.linalg.norm(N7_N2) <= max_bond_distance and np.linalg.norm(O6_N1) <= max_bond_distance):
                                distanceDict.get(keys)["," + key + "," + key2 + ",N7_N2"] = N7_N2
                                distanceDict.get(keys)["," + key + "," + key2 + ",O6_N1"] = O6_N1

                # used for limitation between neighbour guainins, need comma at the begining for this


    return distanceDict




def sort_list(list1, list2):
    zipped_pairs = zip(list2, list1)

    z = [x for _, x in sorted(zipped_pairs)]

    return z




BondAnglesDict= {}

def findTetrads_AngleMethod(angleDict):
    tetrads = {}

    for keys in angleDict.keys():
        tetrads[keys] = []
        for key in angleDict[keys]:

            line = key.split(",")
            check = True
            i = 0

            bond = line[5]
            next_bond = ""
            corresponding_bond = ""
            reverse_corresponding_bond = ""


            # determines the bond that must be established for a specific pair
            if (bond.__contains__("N2_N7")):
                next_bond = "N7_N2"
                corresponding_bond = "N1_O6"
                reverse_corresponding_bond = "O6_N1"
            elif (bond.__contains__("N7_N2")):
                next_bond = "N2_N7"
                corresponding_bond = "O6_N1"
                reverse_corresponding_bond = "N1_O6"
            elif (bond.__contains__("O6_N1")):
                next_bond = "N1_O6"
                corresponding_bond = "N7_N2"
                reverse_corresponding_bond = "N2_N7"
            elif (bond.__contains__("N1_O6")):
                next_bond = "O6_N1"
                corresponding_bond = "N2_N7"
                reverse_corresponding_bond = "N7_N2"

            # determines the guanins in g4 by checking the bonds between them
            for tetrad_Guanines in tetrads[keys]:

                # if there is one of the guanine and not the other, and if there is other corresponding bond and distance of the corresponding bond is less than 5
                if (int(line[1]) in tetrad_Guanines and not (int(line[3]) in tetrad_Guanines)
                        and (any(("," + line[3] + ",DG," + line[1] + ",DG," + next_bond) in x for x in
                                angleDict[keys].keys()) and
                             any(("," + line[1] + ",DG," + line[3] + ",DG," + corresponding_bond) in x for x in
                                 angleDict[keys].keys()) and
                             any(("," + line[3] + ",DG," + line[1] + ",DG," + reverse_corresponding_bond) in x for x in
                                 angleDict[keys].keys())
                        )):
                    check = False
                    tetrads[keys][i].append(int(line[3]))

                elif (not (int(line[1]) in tetrad_Guanines) and int(line[3]) in tetrad_Guanines
                      and (any(("," + line[3] + ",DG," + line[1] + ",DG," + next_bond) in x for x in
                              angleDict[keys].keys()) and
                           any(("," + line[1] + ",DG," + line[3] + ",DG," + corresponding_bond) in x for x in
                               angleDict[keys].keys()) and
                           any(("," + line[3] + ",DG," + line[1] + ",DG," + reverse_corresponding_bond) in x for x in
                               angleDict[keys].keys())
                      )):
                    check = False
                    tetrads[keys][i].append(int(line[1]))

                elif (int(line[1]) in tetrad_Guanines and int(line[3]) in tetrad_Guanines):
                    check = False
                i += 1

            #to prevent appending same guanine again
            for tetrad_Guanines_2 in tetrads[keys]:
                if ((int(line[1]) in tetrad_Guanines_2 or int(line[3]) in tetrad_Guanines_2) and len(
                        tetrad_Guanines_2) == 4):
                    check = False

            # since there is no guanine in tetrad list at first iteration, it will firstly add the first guanine pair two tetrad list here before entering previous loop
            if ( (any(("," + line[3] + ",DG," + line[1] + ",DG," + next_bond) in x for x in angleDict[keys].keys()) and
                                        any(("," + line[1] + ",DG," + line[3] + ",DG," + corresponding_bond) in x for x in angleDict[keys].keys()) and
                                        any(("," + line[3] + ",DG," + line[1] + ",DG," + reverse_corresponding_bond) in x for x in angleDict[keys].keys())
            )  ):# check.__eq__(True) and excluded
                tetrads[keys].append([int(line[1]), int(line[3])])

    return tetrads





def checkTetradCycle(tetrads, angleDict_narrowRange, angleDict_wideRange):
    wideRangeList = []
    for keys in tetrads.keys():

        tetrad_no = -1
        delete_list = []


        for tetrad in tetrads[keys]:
            tetrad_no += 1

            loop = isLoop(tetrad, keys, angleDict_narrowRange, tetrad_no, tetrads)

            if (not loop):
                loop = isLoop(tetrad, keys, angleDict_wideRange, tetrad_no, tetrads)
                if(loop):
                    wideRangeList.append(keys)
                    wideRangeList.append(str(tetrad))

                if(not loop):
                    delete_list.append(tetrad_no)

        j = 0
        for tetradNo in delete_list:
            if (j > 0):
                tetradNo -= j
            del (tetrads[keys][tetradNo])
            j += 1

    return wideRangeList


def isLoop(tetrad, keys, angleOrDistanceDict, tetrad_no, tetrads):
    loop= False
    for g1 in tetrad:
        if (loop):
            break
        path = [0, 0, 0, 0]
        path[0] = int(g1)

        for g2 in tetrad:
            if (loop):
                break

            path[1] = int(g2)

            for g3 in tetrad:
                if (loop):
                    break

                path[2] = int(g3)

                for g4 in tetrad:
                    if (loop):
                        break

                    path[3] = int(g4)

                    if (((
                            (any(("," + str(g1) + ",DG," + str(g2) + ",DG,N2_N7") in x for x in angleOrDistanceDict[keys].keys()) and
                            any(("," + str(g2) + ",DG," + str(g3) + ",DG,N2_N7") in x for x in angleOrDistanceDict[keys].keys()) and
                            any(("," + str(g3) + ",DG," + str(g4) + ",DG,N2_N7") in x for x in angleOrDistanceDict[keys].keys()) and
                            any(("," + str(g4) + ",DG," + str(g1) + ",DG,N2_N7") in x for x in angleOrDistanceDict[keys].keys()))#only one once upon a time
                            and
                            (any(("," + str(g1) + ",DG," + str(g2) + ",DG,N1_O6") in x for x in angleOrDistanceDict[keys].keys()) and
                             any(("," + str(g2) + ",DG," + str(g3) + ",DG,N1_O6") in x for x in angleOrDistanceDict[keys].keys()) and
                             any(("," + str(g3) + ",DG," + str(g4) + ",DG,N1_O6") in x for x in angleOrDistanceDict[keys].keys()) and
                             any(("," + str(g4) + ",DG," + str(g1) + ",DG,N1_O6") in x for x in angleOrDistanceDict[keys].keys()))
                    ) or
                            (any(("," + str(g1) + ",DG," + str(g2) + ",DG,N7_N2") in x for x in angleOrDistanceDict[keys].keys()) and
                            any(("," + str(g2) + ",DG," + str(g3) + ",DG,N7_N2") in x for x in angleOrDistanceDict[keys].keys()) and
                            any(("," + str(g3) + ",DG," + str(g4) + ",DG,N7_N2") in x for x in angleOrDistanceDict[keys].keys()) and
                            any(("," + str(g4) + ",DG," + str(g1) + ",DG,N7_N2") in x for x in angleOrDistanceDict[keys].keys()))
                            and
                            (any(("," + str(g1) + ",DG," + str(g2) + ",DG,O6_N1") in x for x in angleOrDistanceDict[keys].keys()) and
                             any(("," + str(g2) + ",DG," + str(g3) + ",DG,O6_N1") in x for x in angleOrDistanceDict[keys].keys()) and
                             any(("," + str(g3) + ",DG," + str(g4) + ",DG,O6_N1") in x for x in angleOrDistanceDict[keys].keys()) and
                             any(("," + str(g4) + ",DG," + str(g1) + ",DG,O6_N1") in x for x in angleOrDistanceDict[keys].keys()))
                            )


                            and
                            len(set([g1, g2, g3, g4])) == 4):

                        loop = True
                        tetrads[keys][tetrad_no] = [g1, g2, g3, g4]
                        break

    return loop


def deleteImproperTetrads(tetrads):

    # delete the possible tetrad whose length is not 4, because these ones can not form a tetrad
    for key in tetrads.keys():
        i = 0
        files = []
        for tetrad in tetrads[key]:

            if (len(tetrad) != 4):
                files.append(i)
            i += 1

        j = 0
        for a in files:
            if (j > 0):
                a -= j
            del (tetrads[key][a])

            j += 1

    files = []

    # delete duplicate g4s
    for key in tetrads.keys():
        foundG4s = {}
        delete = []
        g4s = tetrads[key]
        index = -1
        for g4 in g4s:
            index += 1
            for element in range(0, len(g4)):

                # if we see a guanine for the first time then add to the dict, if we see the same guanine for the second time than delete this tetrad
                if (not any(("," + str(g4[element]) + ",") in x for x in foundG4s.keys())):
                    foundG4s["," + str(g4[element]) + ","] = "," + str(g4[element]) + ","
                elif (any(("," + str(g4[element]) + ",") in x for x in foundG4s.keys())):
                    delete.append(index)
                    break
        j = 0
        for element in delete:
            if (j > 0):
                element -= j
            del (tetrads[key][element])
            j += 1

    # delete g4s with less than 2 tetrads
    for key in tetrads.keys():
        if (len(tetrads[key]) < 2):
            files.append(key)

    for element in files:
        del (tetrads[element])



# from top to bottom vertically
def orderTetrads(tetrads, atomsDNA):
    # order the tetrads from top to bottom vertically
    ordered_tetrads = deepcopy(tetrads)
    for key in ordered_tetrads.keys():
        order = []  # np.array()
        g4s = ordered_tetrads[key]

        # find the middle points of all tetrads
        for g4 in g4s:

            middle_point = np.array([0.0, 0.0, 0.0])

            # find one tetrad' s middle point
            for point in range(0, len(g4)):

                difference = (np.array(
                    list(map(float, (atomsDNA[key][str(g4[point]) + ",DG"]["N7"])))) + np.array(
                    list(map(float, (atomsDNA[key][str(g4[point]) + ",DG"]["N1"]))))) / 2
                middle_point += difference

            middle_point = middle_point / 4  # count
            order.append(middle_point)

        max = abs(np.linalg.norm(order[0] - order[1]))
        index1 = 0
        index2 = 1
        counter1 = -1

        # find the tetreads which have highest distance between them
        for a in order:
            counter1 += 1
            counter2 = -1
            for b in order:
                counter2 += 1
                distance = abs(np.linalg.norm(a - b))  # find the distance
                if (distance > max):
                    max = distance
                    index1 = counter1
                    index2 = counter2

        count = len(order)
        start = 0
        index = 1
        order2 = []

        order2.append(index1)

        min = abs(np.linalg.norm(order[0] - order[1]))


        # find the closest tetrads
        while (len(order2) != count):
            last_element = order2[len(order2) - 1]
            min = 1000

            for a in range(0, count):
                if (np.array_equal(order[last_element], order[a]) or a in order2):
                    continue
                if (abs(np.linalg.norm(order[a] - order[last_element])) < min):
                    min = abs(np.linalg.norm(order[a] - order[last_element]))
                    index = a

            order2.append(index)

        # use tuple since list is updated and we need initial list
        temp = tuple(ordered_tetrads[key])
        # print(temp)
        counter_ = -1
        for no in order2:
            counter_ += 1
            ordered_tetrads[key][counter_] = temp[no]


    return ordered_tetrads


# numerically
def sortSequenceInLoop(stack, sequence, key, reversed):
    sort = []

    first = int(sequence[len(sequence) - 1])
    second = int(sequence[len(sequence) - 2])
    sort.append(first)
    sort.append(second)
    if (len(stack) == 2):
        sort.sort(reverse=reversed)
        sequence[len(sequence) - 2] = str(sort[0])
        sequence[len(sequence) - 1] = str(sort[1])
    elif (len(stack) == 3):

        third = int(sequence[len(sequence) - 3])
        sort.append(third)
        sort.sort(reverse=reversed)
        sequence[len(sequence) - 3] = str(sort[0])
        sequence[len(sequence) - 2] = str(sort[1])
        sequence[len(sequence) - 1] = str(sort[2])
    elif (len(stack) == 4):
        third = int(sequence[len(sequence) - 3])
        forth = int(sequence[len(sequence) - 4])
        sort.append(third)
        sort.append(forth)
        sort.sort(reverse=reversed)
        sequence[len(sequence) - 4] = str(sort[0])
        sequence[len(sequence) - 3] = str(sort[1])
        sequence[len(sequence) - 2] = str(sort[2])
        sequence[len(sequence) - 1] = str(sort[3])

# find and return loops
def findLoops(tetrads, narrowRangeAngle_or_DistanceDict, wideRangeAngle_or_DistanceDict, wideRangeList, atoms,reversed):
    # find loops
    loops = {}
    if(len(wideRangeList)>0):
        print("wideeee",wideRangeList)

    for key in tetrads:
        for tetrad in tetrads[key]:
            tetrad.sort()

    ordered_tetrads = orderTetrads(tetrads, atoms)


    for key in tetrads.keys():


        angle_or_DistanceDict = narrowRangeAngle_or_DistanceDict
        if (wideRangeList is not None and key in wideRangeList):
            angle_or_DistanceDict = wideRangeAngle_or_DistanceDict


        length = len(tetrads[key])  # shows how many tetrads there are in the structure
        g4s = ordered_tetrads[key]
        foundG4s = {}
        first_index = -1
        second_index = 0
        index = 0
        loops[key] = {}

        # assign each guanin in tetrad to the index it resides in to a map
        for g4 in g4s:
            first_index += 1
            #g4Points = g4.split(",")
            for point in range(0, len(g4)):
                foundG4s[g4[point]] = first_index#int maps to int

        min = 1000
        sequence = []
        stack = []

        order = ""
        counter = 0
        max_guanine = 0
        min_guanine = 10000
        choosenGuanine = 0

        # find the last guanin in all tetrads
        for key4 in foundG4s:
            if (int(key4) > max_guanine):
                max_guanine = int(key4)

        firstG = 1000
        # find the first guanin in all tetrads
        for key4 in foundG4s:
            if (int(key4) < min_guanine):
                min_guanine = int(key4)

        first_guanine = 0
        if(reversed):
            first_guanine = max_guanine  # tetrads[key][0][0]#.split(",")[1]#int
            choosenGuanine = min_guanine
        else:
            first_guanine = min_guanine  # tetrads[key][0][0]#.split(",")[1]#int
            choosenGuanine = max_guanine

        # append first guanin to sequence
        sequence.append(str(first_guanine))
        # append the finst guanin's index to stack
        stack.append(foundG4s[first_guanine])  # list(foundG4s)[0]

        while (True):

            if (str(choosenGuanine) in sequence):
                loops[key] = sequence
                break

            guanin = ""
            min = 1000


            # find the closest guanin to the last member (guanin) of our sequence
            for key2 in foundG4s:

                if (str(key2) in sequence):
                    continue
                difference = abs(int(sequence[len(sequence) - 1]) - int(key2))
                if (difference < min):  # and not stack.__contains__(foundG4s[key3])
                    min = difference
                    guanin = str(key2)

            sequence.append(guanin)
            stack.append(foundG4s[int(guanin)])

            # use tuple since normal list is updated for upcoming calculations, tuple is not affected
            pre = deepcopy(stack)
            sorted_list = sorted(stack)
            reversed_list = sorted(stack, reverse=True)
            stack = pre


            closest_guanin = ""
            enter = ""

            # if we are at last guainin, we may have again a loop
            if (sequence.__contains__(str(choosenGuanine))):

                if (( len(stack) > 2 and (len(stack) != len(set(stack)))) or (len(stack) > 2 and (
                        not (stack == list(dict.fromkeys(sorted_list))
                            or stack == list(dict.fromkeys(reversed_list)))))):

                    sortSequenceInLoop(stack, sequence, key, reversed)
                    lengthStack = len(stack)
                    stack = []
                    decideLoopType(angle_or_DistanceDict, key, sequence, stack, foundG4s, tetrads, lengthStack, length, atoms)
                elif(len(stack) == 2 and stack[0] == stack[1]):
                    sortSequenceInLoop(stack, sequence,key, reversed)
                    stack = []
                    decideLoopTypeForTwoGs(angle_or_DistanceDict, key, sequence, stack, foundG4s, tetrads)
                loops[key] = sequence
                break


            if (len(stack) != len(set(stack))):
                enter = "half"

            for a in range(0, length):

                if (a in stack):
                    enter = "full"
                else:
                    enter = ""
                    break

            # both guanins belong to same index
            if (len(stack) == 2 and stack[0] == stack[1]):


                if (str(choosenGuanine) in sequence):
                    loops[key] = sequence
                    break

                sortSequenceInLoop(stack, sequence,key, reversed)
                stack = []
                decideLoopTypeForTwoGs(angle_or_DistanceDict, key, sequence, stack, foundG4s, tetrads)



            elif (enter.__eq__("full") and (stack == list(dict.fromkeys(sorted_list))
                            or stack == list(dict.fromkeys(reversed_list))) ) :  # after all g4s-        len(stack) == length    enter.__contains__("full")

                sortSequenceInLoop(stack, sequence, key, reversed)  # length can be used instaed of len(stack)

                min = 1000
                closest_guanin = 0

                # find the closest guanin to the last one
                for key3 in foundG4s:

                    if (str(key3) in sequence):
                        continue
                    difference = abs(int(sequence[len(sequence) - 1]) - int(key3))

                    if (difference < min):  # and not stack.__contains__(foundG4s[key3])
                        min = difference
                        closest_guanin = key3

                last_guanin = int(sequence[len(sequence) - 1])
                last_second_guanin = int(sequence[len(sequence) - 2])

                last_guanin_N7 = [float(i) for i in atoms[key][str(last_guanin) + ",DG"]["N7"]]
                last_guanin_N1 = [float(i) for i in atoms[key][str(last_guanin) + ",DG"]["N1"]]

                last_second_guanin_N7 = [float(i) for i in atoms[key][str(last_second_guanin) + ",DG"]["N7"]]
                last_second_guanin_N1 = [float(i) for i in atoms[key][str(last_second_guanin) + ",DG"]["N1"]]

                closest_guanin_N7 = [float(i) for i in atoms[key][str(closest_guanin) + ",DG"]["N7"]]
                closest_guanin__N1 = [float(i) for i in atoms[key][str(closest_guanin) + ",DG"]["N1"]]

                last_guanin_center = (np.array(last_guanin_N7) + np.array(last_guanin_N1)) / 2
                last_second_guanin_center = (np.array(last_second_guanin_N7) + np.array(
                    last_second_guanin_N1)) / 2
                closest_guanin_center = (np.array(closest_guanin_N7) + np.array(
                    closest_guanin__N1)) / 2

                closestGuanines = []


                if (foundG4s[last_guanin] !=foundG4s[last_second_guanin] and
                        (abs(np.linalg.norm(last_guanin_center - closest_guanin_center)) <
                         abs(np.linalg.norm(last_guanin_center - last_second_guanin_center))
                        )):

                    if ((any(("," + str(last_guanin) + ",DG," + str(last_second_guanin) + ",DG,") in x for x in
                            angle_or_DistanceDict[key].keys()) or
                            any(("," + str(last_second_guanin) + ",DG," + str(last_guanin) + ",DG,") in x for x in
                                angle_or_DistanceDict[key].keys())) and foundG4s[last_second_guanin] == foundG4s[last_guanin]):

                        stack = []
                        sequence.insert(len(sequence) - 1, "reversal")
                        sequence.append(str(closest_guanin))
                        stack.append(foundG4s[last_guanin])
                        stack.append(foundG4s[closest_guanin])

                else:

                    if ((any(("," + str(last_guanin) + ",DG," + str(closest_guanin) + ",DG,") in x for x in
                            angle_or_DistanceDict[key].keys()) or
                            any(("," + str(closest_guanin) + ",DG," + str(last_guanin) + ",DG,") in x for x in
                                angle_or_DistanceDict[key].keys())) and foundG4s[closest_guanin] == foundG4s[last_guanin]):
                        stack = []
                        sequence.append("lateral")
                        sequence.append(str(closest_guanin))
                        stack.append(foundG4s[closest_guanin])
                    elif (foundG4s[closest_guanin] == foundG4s[last_guanin]):

                        stack = []
                        sequence.append("diagonal")
                        sequence.append(str(closest_guanin))
                        stack.append(foundG4s[closest_guanin])

                    else:
                        stack = []
                        sequence.append("reversal")  # sequence[len(sequence) - 1] = "reversal"  ## recent indentation
                        sequence.append(str(closest_guanin))
                        stack.append(foundG4s[closest_guanin])


            elif ((len(stack) > 2 and (len(stack) != len(set(stack)))) or (len(stack) > 2 and (not (stack == list(dict.fromkeys(sorted_list))
                            or stack == list(dict.fromkeys(reversed_list)))))):

                if (str(choosenGuanine) in sequence):
                    loops[key] = sequence
                    break


                sortSequenceInLoop(stack, sequence,key, reversed)
                lengthStack = len(stack)
                stack = []
                decideLoopType(angle_or_DistanceDict, key, sequence, stack, foundG4s, tetrads, lengthStack, length, atoms)

    return loops


def  crossProduct(p0, p1, p2):
    x0, x1, x2 = [float(x) for x in p0]
    y0, y1, y2 = [float(x) for x in p1]
    z0, z1, z2 = [float(x) for x in p2]

    A = np.array([x0,x1,x2])
    B = np.array([y0,y1,y2])
    C = np.array([z0,z1,z2])

    result2 = np.dot(getNormal(A,B,C),(np.cross(B-A,C-A)))
    return result2

def ONZ(tetrads, narrowRangeAngle_or_DistanceDict, wideRangeAngle_or_DistanceDict, wideRangeList, atoms):
    print("ONZ Classification")

    for key in tetrads:
        for tetrad in tetrads[key]:
            tetrad.sort()

    ordered_tetrads = orderTetrads(tetrads, atoms)

    ONZ = {}
    tetrad_onz = {}
    for key in tetrads.keys():

        angle_or_DistanceDict = narrowRangeAngle_or_DistanceDict
        if (wideRangeList is not None and key in wideRangeList):
            angle_or_DistanceDict = wideRangeAngle_or_DistanceDict


        for tetrad in tetrads[key]:


            if(not key in ONZ.keys()):
                # print(key,"-")

                ONZ[key] = []
                tetrad_onz[key] = []


            if((any(("," + str(tetrad[0]) + ",DG," + str(tetrad[1]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                  ) and
            any(("," + str(tetrad[1]) + ",DG," + str(tetrad[2]) + ",DG,") in y for y in angle_or_DistanceDict[key].keys()
                ) and
            any(("," + str(tetrad[2]) + ",DG," + str(tetrad[3]) + ",DG,") in z for z in angle_or_DistanceDict[key].keys()
                ) and
            any(("," + str(tetrad[3]) + ",DG," + str(tetrad[0]) + ",DG,") in t for t in angle_or_DistanceDict[key].keys()
                ) )
            or
            (any(("," + str(tetrad[0]) + ",DG," + str(tetrad[3]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                     ) and
            any(("," + str(tetrad[3]) + ",DG," + str(tetrad[2]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                ) and
            any(("," + str(tetrad[2]) + ",DG," + str(tetrad[1]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                ) and
            any(("," + str(tetrad[1]) + ",DG," + str(tetrad[0]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                ) )
            ):
                ONZ[key].append("O")

                tetrad_onz[key].append(tetrad)
            elif((any(("," + str(tetrad[0]) + ",DG," + str(tetrad[1]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                     ) and
            any(("," + str(tetrad[1]) + ",DG," + str(tetrad[3]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                ) and
            any(("," + str(tetrad[3]) + ",DG," + str(tetrad[2]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                ) and
            any(("," + str(tetrad[2]) + ",DG," + str(tetrad[0]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                ) )
            or
            (any(("," + str(tetrad[0]) + ",DG," + str(tetrad[2]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                     ) and
            any(("," + str(tetrad[2]) + ",DG," + str(tetrad[3]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                ) and
            any(("," + str(tetrad[3]) + ",DG," + str(tetrad[1]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                ) and
            any(("," + str(tetrad[1]) + ",DG," + str(tetrad[0]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                ))
            ):
                ONZ[key].append("N")


                tetrad_onz[key].append(tetrad)
            elif((any(("," + str(tetrad[0]) + ",DG," + str(tetrad[3]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                    ) and
            any(("," + str(tetrad[3]) + ",DG," + str(tetrad[1]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                ) and
            any(("," + str(tetrad[1]) + ",DG," + str(tetrad[2]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                ) and
            any(("," + str(tetrad[2]) + ",DG," + str(tetrad[0]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()))
            or
            (any(("," + str(tetrad[0]) + ",DG," + str(tetrad[2]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                     ) and
            any(("," + str(tetrad[2]) + ",DG," + str(tetrad[1]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                ) and
            any(("," + str(tetrad[1]) + ",DG," + str(tetrad[3]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                ) and
            any(("," + str(tetrad[3]) + ",DG," + str(tetrad[0]) + ",DG,") in x for x in angle_or_DistanceDict[key].keys()
                ) )):
                ONZ[key].append("Z")

                tetrad_onz[key].append(tetrad)


        print(ONZ[key])
        print(tetrad_onz[key])
        print("-----------")
    return ONZ, tetrad_onz


def decideLoopTypeForTwoGs(distancOrAngleeDict, key, sequence, stack, foundG4s, tetrads):
    last_guanin = int(sequence[len(sequence) - 1])
    last_second_guanin = int(sequence[len(sequence) - 2])

    # check if these guains are present in previously found pairs
    if ((any(("," + str(last_guanin) + ",DG," + str(last_second_guanin) + ",DG,") in x for x in
            distancOrAngleeDict[key].keys()) or
            any(("," + str(last_second_guanin) + ",DG," + str(last_guanin) + ",DG,") in x for x in
                distancOrAngleeDict[key].keys())) and last_second_guanin in tetrads[key][foundG4s[last_guanin]]):
        temp = sequence[len(sequence) - 1]
        sequence[len(sequence) - 1] = "lateral"
        sequence.append(temp)
        stack.append(foundG4s[int(temp)])
    elif (last_second_guanin in tetrads[key][foundG4s[last_guanin]]):

        temp = sequence[len(sequence) - 1]
        sequence[len(sequence) - 1] = "diagonal"
        sequence.append(temp)
        stack.append(foundG4s[int(temp)])


# reversal, diagona or lateral
def decideLoopType(distancOrAngleeDict, key, sequence, stack, foundG4s, tetrads, lengthStack, lengthTetrads, atoms):
    last_guanin = int(sequence[len(sequence) - 1])
    last_second_guanin = int(sequence[len(sequence) - 2])
    last_third_guanin = int(sequence[len(sequence) - 3])

    last_guanin_center_N7 = [float(i) for i in atoms[key][str(last_guanin)+",DG"]["N7"]]
    last_guanin_center_N1 = [float(i) for i in atoms[key][str(last_guanin)+",DG"]["N1"]]

    last_second_guanin_center_N7 = [float(i) for i in atoms[key][str(last_second_guanin)+",DG"]["N7"]]
    last_second_guanin_center_N1 = [float(i) for i in atoms[key][str(last_second_guanin)+",DG"]["N1"]]

    last_third_guanin_center_N7 = [float(i) for i in atoms[key][str(last_third_guanin)+",DG"]["N7"]]
    last_third_guanin_center_N1 = [float(i) for i in atoms[key][str(last_third_guanin)+",DG"]["N1"]]

    last_guanin_center = (np.array(last_guanin_center_N7) + np.array(last_guanin_center_N1)) / 2
    last_second_guanin_center = (np.array(last_second_guanin_center_N7) + np.array(last_second_guanin_center_N1)) / 2
    last_third_guanin_center = (np.array(last_third_guanin_center_N7) + np.array(last_third_guanin_center_N1)) / 2

    if (foundG4s[last_guanin] == foundG4s[last_second_guanin]):

        if ((any(("," + str(last_guanin) + ",DG," + str(last_second_guanin) + ",DG,") in x for x in
                 distancOrAngleeDict[key].keys()) or
             any(("," + str(last_second_guanin) + ",DG," + str(last_guanin) + ",DG,") in x for x in
                 distancOrAngleeDict[key].keys())) and foundG4s[last_guanin] == foundG4s[last_second_guanin]):

            sequence[len(sequence) - 1] = "lateral"
            sequence.append(str(last_guanin))
            stack.append(foundG4s[last_guanin])
        elif (foundG4s[last_guanin] == foundG4s[last_second_guanin]):

            sequence[len(sequence) - 1] = "diagonal"
            sequence.append(str(last_guanin))
            stack.append(foundG4s[last_guanin])

    elif(foundG4s[last_guanin] != foundG4s[last_second_guanin] and
        (
                    abs(np.linalg.norm(last_second_guanin_center- last_guanin_center)) <
                    abs(np.linalg.norm(last_second_guanin_center- last_third_guanin_center))
        )):
        sequence.insert(len(sequence) - 2, "reversal")
        stack.append(foundG4s[last_second_guanin])
        stack.append(foundG4s[last_guanin])
    else:
        sequence.insert(len(sequence) - 1, "reversal")
        stack.append(foundG4s[last_guanin])


RED = '\033[31m'
ENDC = '\033[m'
GREEN = '\033[32m'

def listG4Structure(pdbName, tetrads, loops, atoms_sequence):
    #print(pdbName)
    print(tetrads[pdbName])
    print(loops[pdbName])
    bgm_list = []
    gfl_list = []
    gf0_list = []
    eight_og_list = []
    gf2_list = []
    lcg_list = []
    zerog_list = []


    enter_bgm = False
    enter_gfl = False
    enter_gf0 = False
    enter_eight_og = False
    enter_gf2 = False
    enter_lcg = False
    enter_zerog = False


    for tetrad in tetrads[pdbName]:

        for index in tetrad:
            if (atoms_sequence[pdbName][str(index)].__contains__("GM")):
                enter_bgm = True
                bgm_list.append(str(index))
            if (atoms_sequence[pdbName][str(index)].__contains__("OG")):
                enter_eight_og = True
                eight_og_list.append(str(index))
            if (atoms_sequence[pdbName][str(index)].__contains__("FL")):
                enter_gfl = True
                gfl_list.append(str(index))
            if (atoms_sequence[pdbName][str(index)].__contains__("F2")):
                enter_gf2 = True
                gf2_list.append(str(index))
            if (atoms_sequence[pdbName][str(index)].__contains__("F0")):
                enter_gf0 = True
                gf0_list.append(str(index))
            if (atoms_sequence[pdbName][str(index)].__contains__("CG")):
                enter_lcg = True
                lcg_list.append(str(index))
            if (atoms_sequence[pdbName][str(index)].__contains__("0G")):
                enter_zerog = True
                zerog_list.append(str(index))

    if (enter_bgm):
        print("BGMs: ", end=" ")
        print(str(bgm_list))
    if (enter_eight_og):
        print("8OGs: ", end=" ")
        print(str(eight_og_list))
    if (enter_gfl):
        print( "GFLs: ", end=" ")
        print( str(gfl_list))
    if (enter_gf2):
        print( "GF2s: ", end=" ")
        print(str(gf2_list))
    if (enter_gf0):
        print( "GF0s: ", end=" ")
        print( str(gf0_list))
    if (enter_lcg):
        print( "LCGs: ", end=" ")
        print(str(lcg_list))
    if (enter_zerog):
        print( "0Gs: ", end=" ")
        print( str(zerog_list))

    # print("loop",loops[pdbName])
    ligand_list = []

    for elements in atoms_sequence[pdbName]:
        if (elements.__contains__("ligand")):
            ligand_list.append(elements)

    for key in atoms_sequence[pdbName]:
        if (atoms_sequence[pdbName][key].__contains__("GM") or atoms_sequence[pdbName][key].__contains__("OG") or atoms_sequence[pdbName][key].__contains__("FL") or
        atoms_sequence[pdbName][key].__contains__("F2") or atoms_sequence[pdbName][key].__contains__("F0") or atoms_sequence[pdbName][key].__contains__("CG")
        or atoms_sequence[pdbName][key].__contains__("0G")):
            if(key+"_ligand" in ligand_list):
                ligand_list.remove(key+"_ligand")

    ligand_enter = False
    if ( len(ligand_list)>0):
        print("Ligands: ", end=" ")
        ligand_enter = True

        print(str(ligand_list))

    return {pdbName: [tetrads[pdbName], [loops[pdbName]]]}, ligand_enter


def bonds(pdbName, distanceDict):
    print(pdbName)
    for element in distanceDict[pdbName]:
        print(element)
        print(distanceDict[pdbName][element])


def calculate(a, b):
    print(abs(np.linalg.norm(a - b)))


def _tetradFromStringToInt(tetradsList):
    tetradListInt = []
    for tetrad in tetradsList:
        element = tetrad.split(",")
        tetradInt = []
        for guanineIndex in range(1, len(element) - 1):
            tetradInt.append(int(element[guanineIndex]))
        tetradInt.sort()
        tetradListInt.append(tetradInt)

    return tetradListInt


def _sortTetradGuanines(tetradsList):
    tetradListInt = []
    for tetrad in tetradsList:
        tetrad.sort()
        tetradListInt.append(tetrad)
    return tetradListInt


def _sortTetradMinToMax(tetradsList):
    for tetrad in range(0, len(tetradsList)):
        minG = 100
        index = tetrad
        copyTetradList = deepcopy(tetradsList)

        for element in range(tetrad, len(tetradsList)):
            if (min(tetradsList[element]) < minG):
                minG = min(tetradsList[element])
                index = element

        tetradsList[tetrad] = copyTetradList[index]
        tetradsList[index] = copyTetradList[tetrad]



def checkTestSet(validationSet, testSet):
    for keysValidation in validationSet:
        if (testSet.keys().__contains__(keysValidation)):
            tetradsValidation = validationSet[keysValidation][0]
            loopValidation = validationSet[keysValidation][1]

            tetradsTest = testSet[keysValidation][0]
            loopTest = testSet[keysValidation][1]

            tetradListIntValidation = _sortTetradGuanines(tetradsValidation)
            _sortTetradMinToMax(tetradsTest)
            _sortTetradMinToMax(tetradListIntValidation)
            print("***")
            print(keysValidation)
            print("Test Set")
            print(tetradsTest)
            print(loopTest)
            print("Validation Set")
            print(tetradListIntValidation)

            for loop in range(0, len(loopValidation)):
                string = loopValidation[loop][0].split(",")
                for element in range(0, len(string)):
                    if (string[element].__contains__("r")):
                        string[element] = "reversal"
                    elif (string[element].__contains__("l")):
                        string[element] = "lateral"
                    elif (string[element].__contains__("d")):
                        string[element] = "diagonal"
                loopValidation[loop] = string

            print(loopValidation)

            print("Tetrad")
            checkTetrad = np.array_equiv(tetradsTest, tetradListIntValidation)
            if (checkTetrad):
                print(checkTetrad)
            else:
                print(RED + str(checkTetrad) + ENDC)

            checkLoop = False
            for element in loopValidation:
                if (np.array_equiv(element, loopTest)):
                    checkLoop = True

            print("Loop")
            if (checkLoop):
                print(checkLoop)
            else:
                print(RED + str(checkLoop) + ENDC)


def convertKeywordsOneComma(from_, to_, data):
    for key in data:
        del_list = []
        add_list = []

        for key2 in data[key]:
            if (key2.__contains__(from_)):
                new_key = key2.replace("," + from_, "," + to_)
                dict = {new_key: data[key][key2]}
                add_list.append(dict)
                del_list.append(key2)
        for key3 in del_list:
            del (data[key][key3])
        for key3 in add_list:
            data[key].update(key3)


def convertKeywordsTwoCommas(from_, to_, data):
    for key in data:
        del_list = []
        add_list = []
        for key2 in data[key]:
            new_key = key2.replace("," + from_ + ",", "," + to_ + ",")
            dict = {new_key: data[key][key2]}
            add_list.append(dict)
            del_list.append(key2)
        for key3 in del_list:
            del (data[key][key3])
        for key3 in add_list:
            data[key].update(key3)


def showDistance(atoms, key, g1, g2, atomOfG1, atomOfG2):
    print(key)
    print(str(g1) + "-" + atomOfG1)
    print(atoms[key][str(g1) + ",DG"][atomOfG1])
    print(str(g2) + "-" + atomOfG2)
    print(atoms[key][str(g2) + ",DG"][atomOfG2])
    a = [float(i) for i in atoms[key][str(g1) + ",DG"][atomOfG1]]
    b = [float(i) for i in atoms[key][str(g2) + ",DG"][atomOfG2]]

    calculate(np.array(a), np.array(b))




def findNormals(pdbList: dict):
    results={}
    for pdbName in pdbList.keys():
        results[pdbName]={}
        nucleotideContent = pdbList[pdbName]
        for nucleotideName in nucleotideContent.keys():
            if nucleotideName.__contains__("G"):
                atomList=nucleotideContent[nucleotideName]
                plane=getNormal(atomList["N1"],atomList["N2"],atomList["N7"])
                results[pdbName][nucleotideName]=plane
    return results


def findHBondGuaninePlaneNormalAngle(hBondVectorList, normalsList, angleTolerance):
    angle_lowerLimit, angle_upperLimit= 90-angleTolerance, 90+angleTolerance
    results = {}

    for pdbName in hBondVectorList.keys():
        results[pdbName] = {}

        for hbondName in hBondVectorList[pdbName].keys():
            hbondVector=  hBondVectorList[pdbName][hbondName]
            donorGname=",".join(hbondName.strip(",").split(",")[:2])
            normalVector= normalsList[pdbName][donorGname]
            angle= findAngleBetween(hbondVector, normalVector) # angle between the first guanine' s normal and bond established by two corresponding guanines

            if(angle>angle_lowerLimit and angle<angle_upperLimit): #60 and 125

                results[pdbName][hbondName]=angle

    return  results


def getNormal(p0,p1,p2):
    # p0, p1, p2 = points
    x0, y0, z0 = [float(x) for x in p0]
    x1, y1, z1 = [float(x) for x in p1]
    x2, y2, z2 = [float(x) for x in p2]

    ux, uy, uz = u = [x1 - x0, y1 - y0, z1 - z0]
    vx, vy, vz = v = [x2 - x0, y2 - y0, z2 - z0]

    u_cross_v = [uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx]

    normal = np.array(u_cross_v)
    return normal

def findAngleBetween(v0,v1):
    unit_vector_1 = v0 / np.linalg.norm(v0)
    unit_vector_2 = v1 / np.linalg.norm(v1)
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    angle = np.arccos(dot_product)
    return np.rad2deg(angle)