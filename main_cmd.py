import sys, CIIS_GQ

# print ('Number of arguments:', len(sys.argv), 'arguments.')
# print ('Argument List:', str(sys.argv))


print()
file = sys.argv[1]
file_name = str(file).split("\\")[-1].split(".")[0]
print("Name of the structure: " + file_name)
atomsDNA, atomsRNA, atomsDNA_sequence, atomsRNA_sequence = CIIS_GQ.readFromFiles([file])
# print(atomsDNA)
# print(atomsRNA)
# print(atomsDNA_sequence)

CIIS_GQ.convertKeywordsOneComma("FL","DG",atomsDNA)
CIIS_GQ.convertKeywordsOneComma("OG","DG",atomsDNA)
CIIS_GQ.convertKeywordsOneComma("F2","DG",atomsDNA)
CIIS_GQ.convertKeywordsOneComma("CG","DG",atomsDNA)
CIIS_GQ.convertKeywordsOneComma("F0","DG",atomsDNA)
CIIS_GQ.convertKeywordsOneComma("GM","DG",atomsDNA)
CIIS_GQ.convertKeywordsOneComma("0G","DG",atomsDNA)

CIIS_GQ.convertKeywordsOneComma("G","DG",atomsDNA)

CIIS_GQ.convertKeywordsOneComma("FL","DG",atomsRNA)
CIIS_GQ.convertKeywordsOneComma("OG","DG",atomsRNA)
CIIS_GQ.convertKeywordsOneComma("F2","DG",atomsRNA)
CIIS_GQ.convertKeywordsOneComma("CG","DG",atomsRNA)
CIIS_GQ.convertKeywordsOneComma("F0","DG",atomsRNA)
CIIS_GQ.convertKeywordsOneComma("GM","DG",atomsRNA)
CIIS_GQ.convertKeywordsOneComma("0G","DG",atomsRNA)

CIIS_GQ.convertKeywordsOneComma("G","DG",atomsRNA)

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

bondVectors = CIIS_GQ.calculateDistances_AngleMethod(atoms, max_bond_distance=5)
guaninePlanes = CIIS_GQ.findNormals(atoms)
angleDict_narrowRange = CIIS_GQ.findHBondGuaninePlaneNormalAngle(bondVectors, guaninePlanes, 40)
angleDict_wideRange = CIIS_GQ.findHBondGuaninePlaneNormalAngle(bondVectors, guaninePlanes, 70)


tetrads_byAngle = {}
tetrads_byAngle = CIIS_GQ.findTetrads_AngleMethod(angleDict_narrowRange)
# print(tetrads_byAngle)
if(any(tetrads_byAngle.values())):
    wideRangeAngleList = CIIS_GQ.checkTetradCycle(tetrads_byAngle, angleDict_narrowRange, angleDict_wideRange)
    CIIS_GQ.deleteImproperTetrads(tetrads_byAngle)
    ordered_tetrads = CIIS_GQ.orderTetrads(tetrads_byAngle, atoms)

    # order tetrads from min to max through from bottom to top
    for key in ordered_tetrads.keys():
        length = len(ordered_tetrads[key])
        if (min(ordered_tetrads[key][length - 1]) < min(ordered_tetrads[key][0])):
            ordered_tetrads[key] = ordered_tetrads[key][::-1]

    loops_byAngle = {}
    loops_byAngle = CIIS_GQ.findLoops(tetrads_byAngle, angleDict_narrowRange, angleDict_wideRange, wideRangeAngleList, atoms, reversed=False)
    ONZ,tetrad_onz = CIIS_GQ.ONZ(ordered_tetrads, angleDict_narrowRange, angleDict_wideRange, wideRangeAngleList, atoms)

    print("Tetrads and Loops: ")
    test_set = {}

    for key in loops_byAngle:
        ligand_enter = False

        dict, ligand_enter = CIIS_GQ.listG4Structure(key, ordered_tetrads, loops_byAngle, atoms_sequence)
        test_set.update(dict)

    CIIS_GQ.drawStructure(test_set, file_name, ONZ, tetrad_onz, atoms_sequence)
else:
    print("It is not a G4 structure")
