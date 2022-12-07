set(AFMODEL_SRC afmodel.c)
set(GEN3DDESC_SRC gen3Ddesc.c)
set(FFGEN_SRC ffgen.c)
set(METDBGEN_SRC metdbgen.c)
set(MFIELDSGEN_SRC mfieldgen.c)
set(GENMOL2ADJMX_SRC gen_mol2_adjmx.c)
set(GENMOLEPOT_SRC molepot.c)
set(GENMOLINTEPOT_SRC molintepot.c)
set(COMPARELEPOT_SRC compare_epot.c)
set(PROTLIGGENMOLEPOT_SRC prot_ligand_mol_epot.c)
set(GENMOLUFFVDWPOT_SRC moluffvdwpot.c)
set(GENMOLENEGPOT_SRC molenegpot.c)
set(GENMOLGPOT_SRC molgpot.c)
set(GENMOLLJPOT_SRC molljpot.c)
set(MOLEPOT3DDESC_SRC molepot3Ddesc.c)
set(MOLEPOT3DSPHERICALDESC_SRC molepot3DSphericalDesc.c)
set(MOLEPOTVS_SRC molepotVS.c)
set(MOLGPOT3DDESC_SRC molgpot3Ddesc.c)
set(MOLGPOT3DSPHERICALDESC_SRC molgpot3DSphericalDesc.c)
set(MOLLJPOT3DDESC_SRC molljpot3Ddesc.c)
set(MOLLJPOT3DSPHERICALDESC_SRC molljpot3DSphericalDesc.c)
set(MAKEVOXELDATABASE_SRC makeVoxelDatabase.c)
set(MOL2ATFP_SRC mol2atfp.c)
set(MOL23DATFP_SRC mol23Datypefp.c)
set(MOL23DEPOTFP_SRC mol23Depotfp.c)
set(MOL23DTWOAFP_SRC mol23Dtwoatomfp.c)
set(MOL23DTWOAOneHotEncoding_SRC mol23dtwoatomOneHotEncoding.c)
set(MOL23DTWOADISTANGLETORSION_SRC mol2twoatomdistangletorsion.c)
set(MOL23DTWOAWEPFP_SRC mol23Dtwoatomwepotfp.c)
set(MOL2DESCGEN_SRC mol2GenDescriptors.c)
set(MOL23DTWOATWFP_SRC mol23Dtwoatomtypewpotfp.c)

set(ALIGN3DONVDWSHAPES_SRC align_3d_on_vdw_shapes.c)
set(ALIGN3DCONFORMERS_SRC align_3d_conformers.c)
set(ALIGN3DPHARMACOPHORE_SRC align_3d_pharmacophore.c)


add_executable(metdbgen ${METDBGEN_SRC})
target_link_libraries(metdbgen libmoldesc)

add_executable(afmodel ${AFMODEL_SRC})
target_link_libraries(afmodel libmoldesc)

add_executable(gen3Ddesc ${GEN3DDESC_SRC})
target_link_libraries(gen3Ddesc libmoldesc)

add_executable(ForceFieldGen ${FFGEN_SRC})
target_link_libraries(ForceFieldGen libmoldesc)

add_executable(MIFieldsGen ${MFIELDSGEN_SRC})
target_link_libraries(MIFieldsGen libmoldesc)

add_executable(GenMol2Adjmx ${GENMOL2ADJMX_SRC})
target_link_libraries(GenMol2Adjmx libmoldesc)

add_executable(GenMolEPot ${GENMOLEPOT_SRC})
target_link_libraries(GenMolEPot libmoldesc)

add_executable(GenMolIntEPot ${GENMOLINTEPOT_SRC})
target_link_libraries(GenMolIntEPot libmoldesc)

add_executable(CmpEPot ${COMPARELEPOT_SRC})
target_link_libraries(CmpEPot libmoldesc)

add_executable(GenProteinLigandMolEPot ${PROTLIGGENMOLEPOT_SRC})
target_link_libraries(GenProteinLigandMolEPot libmoldesc)

add_executable(GenMolUFFVdWPot ${GENMOLUFFVDWPOT_SRC})
target_link_libraries(GenMolUFFVdWPot libmoldesc)

add_executable(GenMolEnegPot ${GENMOLENEGPOT_SRC})
target_link_libraries(GenMolEnegPot libmoldesc)

add_executable(GenMolGPot ${GENMOLGPOT_SRC})
target_link_libraries(GenMolGPot libmoldesc)

add_executable(GenMolLJPot ${GENMOLLJPOT_SRC})
target_link_libraries(GenMolLJPot libmoldesc)

add_executable(MolEPot3DDesc ${MOLEPOT3DDESC_SRC})
target_link_libraries(MolEPot3DDesc libmoldesc)

add_executable(MolGPot3DDesc ${MOLGPOT3DDESC_SRC})
target_link_libraries(MolGPot3DDesc libmoldesc)

add_executable(MolLJPot3DDesc ${MOLLJPOT3DDESC_SRC})
target_link_libraries(MolLJPot3DDesc libmoldesc)

add_executable(MolEPot3DSphericalDesc ${MOLEPOT3DSPHERICALDESC_SRC})
target_link_libraries(MolEPot3DSphericalDesc libmoldesc)

add_executable(MOLEpotVS ${MOLEPOTVS_SRC})
target_link_libraries(MOLEpotVS libmoldesc)

add_executable(MolGPot3DSphericalDesc ${MOLGPOT3DSPHERICALDESC_SRC})
target_link_libraries(MolGPot3DSphericalDesc libmoldesc)

add_executable(MolLJPot3DSphericalDesc ${MOLLJPOT3DSPHERICALDESC_SRC})
target_link_libraries(MolLJPot3DSphericalDesc libmoldesc)

add_executable(makeMolVoxelDB ${MAKEVOXELDATABASE_SRC})
target_link_libraries(makeMolVoxelDB libmoldesc)

add_executable(Mol2AtomTypeFP ${MOL2ATFP_SRC})
target_link_libraries(Mol2AtomTypeFP libmoldesc)

add_executable(Mol23DAtomTypeFP ${MOL23DATFP_SRC})
target_link_libraries(Mol23DAtomTypeFP libmoldesc)

add_executable(Mol23DAtomEpotFP ${MOL23DEPOTFP_SRC})
target_link_libraries(Mol23DAtomEpotFP libmoldesc)

add_executable(Mol23DTwoAtomFP ${MOL23DTWOAFP_SRC})
target_link_libraries(Mol23DTwoAtomFP libmoldesc)

add_executable(Mol23DTwoAtomOneHotEncoding ${MOL23DTWOAOneHotEncoding_SRC})
target_link_libraries(Mol23DTwoAtomOneHotEncoding libmoldesc)

add_executable(Mol23DTwoAtomDistanceAngleTorsion ${MOL23DTWOADISTANGLETORSION_SRC})
target_link_libraries(Mol23DTwoAtomDistanceAngleTorsion libmoldesc)

add_executable(Mol23DTwoAtomEpotFP ${MOL23DTWOAWEPFP_SRC})
target_link_libraries(Mol23DTwoAtomEpotFP libmoldesc)

add_executable(Mol2DescriptorGenerator ${MOL2DESCGEN_SRC})
target_link_libraries(Mol2DescriptorGenerator libmoldesc)

add_executable(Mol23DTwoAtomWFP ${MOL23DTWOATWFP_SRC})
target_link_libraries(Mol23DTwoAtomWFP libmoldesc)

add_executable(Align3DMolsOnVDWShapes ${ALIGN3DONVDWSHAPES_SRC})
target_link_libraries(Align3DMolsOnVDWShapes libmoldesc)

add_executable(Align3DConformers ${ALIGN3DCONFORMERS_SRC})
target_link_libraries(Align3DConformers libmoldesc)

add_executable(Align3DPharmacophore ${ALIGN3DPHARMACOPHORE_SRC})
target_link_libraries(Align3DPharmacophore libmoldesc)
