#-----------------------------------------
#          SlicerCervicalSpine master
#          updated: 16.6.2019
#-----------------------------------------
cmake_minimum_required(VERSION 3.5)

project(SlicerCervicalSpine)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/SlicerCervicalSpine")
set(EXTENSION_CATEGORY "VisSimTools")
set(EXTENSION_CONTRIBUTORS "Ibraheem Al-Dhamari (VisSim, Koblenz-Landau University)")
set(EXTENSION_DESCRIPTION "The Cervical spine is the upper part of the spine that supports the head and its movements. It containes seven vertebrae C1 to C7. There are two modules in this extension (more are coming). Both modules do the same tasks but the first one for a single vertebra and the second module for the complete cervical spine. 

Features:
 1. Work on multi-modal images (so far it is tested against CT and MRI).
 2. fast (around 5 seconds). 2. Automatic segmentation.
 3. Automatic ligaments points detection. 
 4. Can be generalized to other parts of the spine or body. Cervical spien sample datasets are available in Slicer Data Store module. 
You are welcome to contribute by correcting, enhancing, optimising the code or sharing your cervical spine dataset. For questions or problems using this extension please use gitHub issue or post in 3D Slicer forum. For more details please cite the related publications: 
  1. Ibraheem AL-Dhamari, Sabine Bauer, Eva Keller and Dietrich Paulus, Automatic Detection of Cervical Spine Ligaments Origin and Insertion Points. Accepted in The IEEE International Symposium on Biomedical Imaging (ISBI)  2019, Venice, Italy.
  2. Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, (2018), Automatic Multi-modal Cervical Spine Image Atlas Segmentation Using Adaptive Stochastic Gradient Descent, Bildverarbeitung für die Medizin 2018 pp 303-308. ")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCervicalSpine/master/SlicerCervicalSpine.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCervicalSpine/master/Screenshots/CS00.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCervicalSpine/master/Screenshots/CVT01.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCervicalSpine/master/Screenshots/CVT02.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCervicalSpine/master/Screenshots/CVT03.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCervicalSpine/master/Screenshots/CVT04.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCervicalSpine/master/Screenshots/CVT05.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCervicalSpine/master/Screenshots/CST01.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCervicalSpine/master/Screenshots/CST02.png https://raw.githubusercontent.com/MedicalImageAnalysisTutorials/SlicerCervicalSpine/master/Screenshots/CST03.png")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

set(EXTENSION_DEPENDS SlicerCochlea)

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(CervicalVertebraTools)
add_subdirectory(CervicalSpineTools)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
