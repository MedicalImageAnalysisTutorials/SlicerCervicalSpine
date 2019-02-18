# run:  ~/sw/Slicer-4.10.0/Slicer --no-main-window --python-script  SlicerCervicalSpineToolsTest.py
#TODO when called from test, no display
#     remove dialogs

import SlicerCervicalSpineTools

CT = SlicerCervicalSpineTools.SlicerCervicalSpineToolsTest()
CT.testSlicerCervicalSpineTools()
