# run:  ~/sw/Slicer-4.10.0/Slicer --no-main-window --python-script  SlicerCervicalVertebraToolsTest.py
#TODO when called from test, no display
#     remove dialogs

import SlicerCervicalVertebraTools

CT = SlicerCervicalVertebraTools.SlicerCervicalVertebraToolsTest()
CT.testSlicerCervicalVertebraTools()
