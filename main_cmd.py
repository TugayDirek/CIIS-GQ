import sys, GIIS_GQ

# print ('Number of arguments:', len(sys.argv), 'arguments.')
# print ('Argument List:', str(sys.argv))


print()
file = sys.argv[1]
file_name = str(file).split("\\")[-1].split(".")[0]
print("Name of the structure: " + file_name)
atomsDNA, atomsRNA, atomsDNA_sequence, atomsRNA_sequence = g4_depictor.readFromFiles([file])
# print(atomsDNA)
# print(atomsRNA)
# print(atomsDNA_sequence)

g4_depictor.convertKeywordsOneComma("FL","DG",atomsDNA)
g4_depictor.convertKeywordsOneComma("OG","DG",atomsDNA)
g4_depictor.convertKeywordsOneComma("F2","DG",atomsDNA)
g4_depictor.convertKeywordsOneComma("CG","DG",atomsDNA)
g4_depictor.convertKeywordsOneComma("F0","DG",atomsDNA)
g4_depictor.convertKeywordsOneComma("GM","DG",atomsDNA)
g4_depictor.convertKeywordsOneComma("0G","DG",atomsDNA)

g4_depictor.convertKeywordsOneComma("G","DG",atomsDNA)

g4_depictor.convertKeywordsOneComma("FL","DG",atomsRNA)
g4_depictor.convertKeywordsOneComma("OG","DG",atomsRNA)
g4_depictor.convertKeywordsOneComma("F2","DG",atomsRNA)
g4_depictor.convertKeywordsOneComma("CG","DG",atomsRNA)
g4_depictor.convertKeywordsOneComma("F0","DG",atomsRNA)
g4_depictor.convertKeywordsOneComma("GM","DG",atomsRNA)
g4_depictor.convertKeywordsOneComma("0G","DG",atomsRNA)

g4_depictor.convertKeywordsOneComma("G","DG",atomsRNA)

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
tetrads_byAngle = g4_depictor.findTetrads_AngleMethod(angleDict_narrowRange)
# print(tetrads_byAngle)
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

    for key in loops_byAngle:
        ligand_enter = False

        dict, ligand_enter = g4_depictor.listG4Structure(key, ordered_tetrads, loops_byAngle, atoms_sequence)
        test_set.update(dict)

    g4_depictor.drawStructure(test_set, file_name, ONZ, tetrad_onz, atoms_sequence)
else:
    print("It is not a G4 structure")
