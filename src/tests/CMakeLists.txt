set(TestGetAngle_SRC testGetAngle.c)
set(TestDihedralAngle_SRC testDihedral.c)
set(TestReadMol2_SRC testReadMol2.c)
set(TestCalcMolCharge_SRC testCalcMolCharge.c)
set(TestMol2TranslateToNewCentre_SRC testMol2TranslateToNewCentre.c)
set(TestMol2ToPQR_SRC testMol2ToPQR.c)
set(TestGenAdjMX_SRC testGenAdiacenceMatrix.c)
set(TestMol2AtomTypeHash_SRC testMol2AtomTypeHash.c)
set(TestAtomDepth_SRC testGetAtomDepth.c)
set(TestNewSDFMolecule_SRC testNewSDFMolecule.c)
set(TestGetPlanarity_SRC testGetPlanarity.c)
set(TestDescriptors_SRC testDescriptors.c)
set(TestGetShapePoints_SRC testGetShapePoints.c)
set(TestGetEMWMWfromMolecularFormula_SRC  testGetEMWMWfromMolecularFormula.c)
set(TestAlign3DConformations_SRC testAlign3DConformations.c)
set(TestAlign3DPharmacophore_SRC testAlign3DPharmacophore.c)
set(TestAlign3DShapes_SRC testAlign3DShapes.c)
set(TestConfRMSD_SRC testRMSDConformations.c)
set(TestMol2RandRotation_SRC testMol2RandRotation.c)
set(TestSphericalInteractionField_SRC testSphericalInteractionField.c)
set(TestVoxel3DFields_SRC testVoxel3DFields testVoxel3DFields.c)
set(TestGetMol2BoxSize_SRC testGetMol2BoxSize.c)
set(TestGetAllDihedralAngleFromId_SRC testGetAllDihedralFromId.c)
set(TestGetAllDihedralAngle_SRC testGetAllDihedral.c)
set(TestGetAllAngles_SRC testGetAllAngles.c)

if(WIN32)
  message(STATUS "Set mingw options")
  set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS_RELEASE)
  set(CMAKE_EXECUTABLE_SUFFIX ".exe")
  set(CMAKE_STATIC_LIBRARY_SUFFIX ".dll")
endif()

add_executable(testAtomDepth ${TestAtomDepth_SRC})
target_link_libraries(testAtomDepth libmoldesc)
add_test(NAME testAtomDepth COMMAND testAtomDepth ../../../src/tests/Data/Cortisol.mol2 10)

add_executable(testGetAngle ${TestGetAngle_SRC})
target_link_libraries(testGetAngle libmoldesc)
add_test(NAME testGetAngle COMMAND testGetAngle)

add_executable(testDihedral ${TestDihedralAngle_SRC})
target_link_libraries(testDihedral libmoldesc)
add_test(NAME testDihedral COMMAND testDihedral)

add_executable(testReadMol2 ${TestReadMol2_SRC})
target_link_libraries(testReadMol2 libmoldesc)
add_test(NAME testReadMol2 COMMAND testReadMol2 ../../../src/tests/Data/Cortisol.mol2)

add_executable(testCalcMolCharge ${TestCalcMolCharge_SRC})
target_link_libraries(testCalcMolCharge libmoldesc)
add_test(NAME testCalcMolCharge COMMAND testCalcMolCharge ../../../src/tests/Data/Cortisol.mol2)

add_executable(testMol2TranslateToNewCentre ${TestMol2TranslateToNewCentre_SRC})
target_link_libraries(testMol2TranslateToNewCentre libmoldesc)
#add_test(NAME testMol2TranslateToNewCentre COMMAND testMol2TranslateToNewCentre ../../../src/tests/Data/Cortisol.mol2)

add_executable(testMol2ToPQR ${TestMol2ToPQR_SRC})
target_link_libraries(testMol2ToPQR libmoldesc)
#add_test(NAME testMol2ToPQR COMMAND testMol2ToPQR ../../../src/tests/Data/Cortisol.mol2)

add_executable(testGenADJmatrix ${TestGenAdjMX_SRC})
target_link_libraries(testGenADJmatrix libmoldesc)
add_test(NAME testGenADJmatrix COMMAND testGenADJmatrix ../../../src/tests/Data/Cortisol.mol2)

add_executable(testMol2AtomTypeHash ${TestMol2AtomTypeHash_SRC})
target_link_libraries(testMol2AtomTypeHash libmoldesc)
add_test(NAME testMol2AtomTypeHash COMMAND testMol2AtomTypeHash ../../../src/tests/Data/Cortisol.mol2)

add_executable(testNewSDFMolecule ${TestNewSDFMolecule_SRC})
target_link_libraries(testNewSDFMolecule libmoldesc)

add_test(NAME testNewSDFMolecule COMMAND testNewSDFMolecule ../../../src/tests/Data/Imatinib.sdf)

add_executable(testGetPlanarity ${TestGetPlanarity_SRC})
target_link_libraries(testGetPlanarity libmoldesc)
add_test(NAME testGetPlanarity COMMAND testGetPlanarity ../../../src/tests/Data/benzene.mol2)

add_executable(testDescriptors ${TestDescriptors_SRC})
target_link_libraries(testDescriptors libmoldesc)
add_test(NAME testDescriptors COMMAND testDescriptors ../../../src/tests/Data/RT212.mol2)

add_executable(testGetShapePoints ${TestGetShapePoints_SRC})
target_link_libraries(testGetShapePoints libmoldesc)
add_test(NAME testGetShapePoints COMMAND testGetShapePoints ../../../src/tests/Data/2ai9.mol2 2ai9.shapepnt.mol2)

add_executable(testGetEMWMWfromMolecularFormula ${TestGetEMWMWfromMolecularFormula_SRC})
target_link_libraries(testGetEMWMWfromMolecularFormula libmoldesc)
add_test(NAME testGetEMWMWfromMolecularFormula COMMAND testGetEMWMWfromMolecularFormula C12H22O11 0)
#add_test(NAME testGetEMWMWfromMolecularFormula COMMAND testGetEMWMWfromMolecularFormula C12H22O11 1)

add_executable(testAlign3DConformations ${TestAlign3DConformations_SRC})
target_link_libraries(testAlign3DConformations libmoldesc)
add_test(NAME testAlign3DConformations COMMAND testAlign3DConformations ../../../src/tests/Data/Conf002_alanine.mol2  ../../../src/tests/Data/Conf100_alanine.mol2)

add_executable(testAlign3DPharmacophore ${TestAlign3DPharmacophore_SRC})
target_link_libraries(testAlign3DPharmacophore libmoldesc)
add_test(NAME testAlign3DPharmacophore COMMAND testAlign3DPharmacophore ../../../src/tests/Data/A.mol2 ../../../src/tests/Data/B.mol2 ../../../src/tests/Data/A_B.txt)

add_executable(testAlign3DShapes ${TestAlign3DShapes_SRC})
target_link_libraries(testAlign3DShapes libmoldesc)
add_test(NAME testAlign3DShapes COMMAND testAlign3DShapes ../../../src/tests/Data/A.mol2 ../../../src/tests/Data/B.mol2)

add_executable(testRMSDConformations ${TestConfRMSD_SRC})
target_link_libraries(testRMSDConformations libmoldesc)
add_test(NAME testRMSDConformations COMMAND testRMSDConformations ../../../src/tests/Data/Conf002_alanine.mol2 ../../../src/tests/Data/Conf100_alanine.mol2)

add_executable(testMol2RandRotation ${TestMol2RandRotation_SRC})
target_link_libraries(testMol2RandRotation libmoldesc)

add_executable(testSphericalInteractionField ${TestSphericalInteractionField_SRC})
target_link_libraries(testSphericalInteractionField libmoldesc)
#add_test(NAME testSphericalInteractionField COMMAND testInteractionField ../../../src/tests/Data/Conf002_alanine.mol2  ../../../src/tests/Data/Conf100_alanine.mol2)

add_executable(testVoxel3DFields ${TestVoxel3DFields_SRC})
target_link_libraries(testVoxel3DFields libmoldesc)

add_executable(testGetMol2BoxSize ${TestGetMol2BoxSize_SRC})
target_link_libraries(testGetMol2BoxSize libmoldesc)

add_executable(testGetAllDihedralAngleFromId ${TestGetAllDihedralAngleFromId_SRC})
target_link_libraries(testGetAllDihedralAngleFromId libmoldesc)

add_executable(testGetAllDihedralAngle ${TestGetAllDihedralAngle_SRC})
target_link_libraries(testGetAllDihedralAngle libmoldesc)

add_executable(testGetAllAngles ${TestGetAllAngles_SRC})
target_link_libraries(testGetAllAngles libmoldesc)

#add_test(NAME testInteractionField COMMAND testInteractionField ../../../src/tests/Data/Conf002_alanine.mol2  ../../../src/tests/Data/Conf100_alanine.mol2)
