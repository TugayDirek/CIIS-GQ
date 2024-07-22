import sys, g4_depictor

# print ('Number of arguments:', len(sys.argv), 'arguments.')
# print ('Argument List:', str(sys.argv))


print()
def main(pdbName):
    file = sys.argv[1]
    #file = "C:\\Users\\tdirek\\PycharmProjects\\tez\\pdb files\\"+pdbName+".ent"
    file_name = str(file).split("\\")[-1].split(".")[0]
    print("Name of the structure: " + file_name)
    atomsDNA, atomsRNA, atomsDNA_sequence, atomsRNA_sequence = g4_depictor.readFromFiles([file])
    #print("atoms DNA", atomsDNA)
    # print("atoms RNA", atomsRNA)
    # print(atomsDNA_sequence)

    g4_depictor.convertKeywordsOneComma("FL","DG",atomsDNA)
    g4_depictor.convertKeywordsOneComma("OG","DG",atomsDNA)
    g4_depictor.convertKeywordsOneComma("F2","DG",atomsDNA)
    g4_depictor.convertKeywordsOneComma("CG","DG",atomsDNA)
    g4_depictor.convertKeywordsOneComma("F0","DG",atomsDNA)
    g4_depictor.convertKeywordsOneComma("GM","DG",atomsDNA)
    g4_depictor.convertKeywordsOneComma("0G","DG",atomsDNA)
    #g4_depictor.convertKeywordsOneComma("ZO","DG",atomsDNA)

    g4_depictor.convertKeywordsOneComma("G","DG",atomsDNA)

    g4_depictor.convertKeywordsOneComma("FL","DG",atomsRNA)
    g4_depictor.convertKeywordsOneComma("OG","DG",atomsRNA)
    g4_depictor.convertKeywordsOneComma("F2","DG",atomsRNA)
    g4_depictor.convertKeywordsOneComma("CG","DG",atomsRNA)
    g4_depictor.convertKeywordsOneComma("F0","DG",atomsRNA)
    g4_depictor.convertKeywordsOneComma("GM","DG",atomsRNA)
    g4_depictor.convertKeywordsOneComma("0G","DG",atomsRNA)
    #g4_depictor.convertKeywordsOneComma("ZO","DG",atomsRNA)

    g4_depictor.convertKeywordsOneComma("G","DG",atomsRNA)
    #print("1111111",atomsDNA)

    atoms = {}
    atoms_sequence = {}
    if len(atomsDNA.keys())!=0:
        atoms = atomsDNA
        atoms_sequence = atomsDNA_sequence
        print("There is 1 DNA structure.")
    elif len(atomsRNA.keys()) != 0:
        atoms = atomsRNA
        atoms_sequence = atomsRNA_sequence
        print("There is 1 RNA structure.")
    print("-----------")
    angleDict = {}



    bondVectors = g4_depictor.calculateDistances_AngleMethod(atoms, max_bond_distance=5)
    guaninePlanes = g4_depictor.findNormals(atoms)
    angleDict_narrowRange = g4_depictor.findHBondGuaninePlaneNormalAngle(bondVectors, guaninePlanes, 40)
    angleDict_wideRange = g4_depictor.findHBondGuaninePlaneNormalAngle(bondVectors, guaninePlanes, 70)


    tetrads_byAngle = {}
    tetrads_byAngle, tetrads_bonds = g4_depictor.findTetrads_AngleMethod(angleDict_narrowRange)
    #print(tetrads_bonds, "---------")

    tetrads_bonds_dict = {}
    for key in tetrads_bonds.keys():
        for element in tetrads_bonds[key]:
            #tetrads_bonds[key] =
            for element2 in element:
                tetrads_bonds_dict[element2[1]+","+element2[5]] = element2[5]
                tetrads_bonds_dict[element2[3]+","+element2[5]] = element2[5]


    #print(tetrads_bonds_dict, "---------")
    #print(tetrads_byAngle, "---------")

    #print(any(tetrads_byAngle.values()))
    if(any(tetrads_byAngle.values())):
        wideRangeAngleList = g4_depictor.checkTetradCycle(tetrads_byAngle, angleDict_narrowRange, angleDict_wideRange)
        g4_depictor.deleteImproperTetrads(tetrads_byAngle)
        ordered_tetrads = g4_depictor.orderTetrads(tetrads_byAngle, atoms)

        # order tetrads from min to max through from bottom to top
        for key in ordered_tetrads.keys():
            length = len(ordered_tetrads[key])
            if (min(ordered_tetrads[key][length - 1]) < min(ordered_tetrads[key][0])):
                ordered_tetrads[key] = ordered_tetrads[key][::-1]

        loops_byAngle = {}
        loops_byAngle = g4_depictor.findLoops(tetrads_byAngle, angleDict_narrowRange, angleDict_wideRange, wideRangeAngleList, atoms, reversed=False)
        ONZ,tetrad_onz = g4_depictor.ONZ(ordered_tetrads, angleDict_narrowRange, angleDict_wideRange, wideRangeAngleList, atoms)

        print("Tetrads and Loops: ")
        test_set = {}
        #print(atoms)
        #print(loops_byAngle)
        for key in loops_byAngle:
            ligand_enter = False

            dict, ligand_enter = g4_depictor.listG4Structure(key, ordered_tetrads, loops_byAngle, atoms_sequence)
            test_set.update(dict)

            # print(key)
            check_heptad = g4_depictor.addHeptads_2(file_name, atoms, 3.5, ordered_tetrads)  # 3.3 was detected as threshold
            n_of_heptad = len(check_heptad[key].keys())
            #print("///", n_of_heptad, check_heptad)
            if n_of_heptad >1:
                print("There are non G-tetrad in this structure")
                #print(check_heptad)


        # for b in check_heptad[key]:
        #     print("-", b)

        g4_depictor.drawStructure(test_set, file_name, ONZ, tetrad_onz, atoms_sequence,check_heptad)


    else:
        print("It is not unimolecular structure")
        #exit(0)

import glob
import os

if __name__=="__main__":

    file_name = ""
    folder_path = "./pdb files"
    files = []
    # save file addresses to a list
    for filename in glob.glob(os.path.join(folder_path, '*.ent')):
        pdbName = filename.split("\\")[1]
        files.append(filename)

    # added for download
    # file_names = []
    # for filename in files:
    #     pdbName = filename.split("\\")[-1].split(".")[0]  # [1]
    #     pdbName = pdbName[3:7]
    #     file_names.append(pdbName)
    # print("There are " + str(len(files)) + " total files at the beginning.")
    # print(file_names)
    #for pdb_name in file_names:
    main("pdb"+"143d")
    print("//////")
    # exit(0)

    #main("pdb1bub")