#======================================================================================
#  3D Slicer [1] plugin that uses elastix toolbox [2] Plugin for Automatic Cervical   # 
#  SPine Image Segmentation and other useful features extraction[3]                   #
#  More info can be found at [3,4 and 5]                                              #
#  Sample vertebra datasets can be downloaded using Slicer Datastore module           #
#                                                                                     #
#  Contributers: Ibraheem Al-Dhamari,  idhamari@uni-koblenz.de                        #
#  [1] https://www.slicer.org                                                         #
#  [2] http://elastix.isi.uu.nl                                                       #
#  [3] Ibraheem AL-Dhamari, Sabine Bauer, Eva Keller and Dietrich Paulus, (2019),     #
#      Automatic Detection Of Cervical Spine Ligaments Origin And Insertion Points,   #
#      IEEE International Symposium on Biomedical Imaging (ISBI), Venice, Italy.      #
#  [4] Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, (2018), Automatic          #
#      Multi-modal Cervical Spine Image Atlas Segmentation Using Adaptive Stochastic  #
#      Gradient Descent, Bildverarbeitung feur die Medizin 2018 pp 303-308.           #  
#  [5] https://mtixnat.uni-koblenz.de                                                 #
#                                                                                     #
#  Updated: 14.2.2019                                                                 #    
#                                                                                     #  
#======================================================================================

import os, re , datetime, time ,shutil, unittest, logging, zipfile, urllib2, stat,  inspect
import sitkUtils, sys ,math, platform, subprocess  
import numpy as np, SimpleITK as sitk
import vtkSegmentationCorePython as vtkSegmentationCore
from __main__ import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *   
from copy import deepcopy
from collections import defaultdict
from os.path import expanduser
from os.path import isfile
from os.path import basename
from PythonQt import BoolResult
from shutil import copyfile
from decimal import *

import Elastix
import SegmentStatistics
#TODO:

# 1. test user inputs
# 2. test windows
# 3. cleaning 
# 4. remove temps nodes and files 
# 5. test again


# Documentation
# Video tutorials
# Uploas all
# VisSimTools and sample data to google drive
 
#   Registration download all stuff for both registration and segmentation.
#   use registration module and commong functions:  
#   Remove inverted transform, it was using in case switching moving and fixed 
# Later:
# - Checking if all above are needed 
# - Cleaning, optimizing, commenting.  
# - Testing in both Windows and Linux. 
# - Supporting DICOM. 
# - Supporting illegal filename.  
# - Using  SlierElastix binaries.   
# - Visualizing the interimediate steps. 
# 
#  
# Terminology
#  img         : ITK image 
#  imgNode     : Slicer Node
#  imgPath     : wholePath + Filename
#  imgFnm      : Filename without the path and the extension
#  imgFileName : Filename without the path

#===================================================================
#                           Main Class
#===================================================================

class CervicalVertebraTools(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        parent.title = "Cervical Vertebra Tools"
        parent.categories = ["VisSimTools"]
        parent.dependencies = []
        parent.contributors = ["Ibraheem Al-Dhamari" ]
        self.parent.helpText += self.getDefaultModuleDocumentationLink()
        #TODO: add sponsor
        parent.acknowledgementText = " This work is sponsored by ................ "
        self.parent = parent
  #end def init
#end class vertebraSeg

    
#===================================================================
#                           Main Widget
#===================================================================
class CervicalVertebraToolsWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    print(" ")
    print("=======================================================")   
    print("   VisSim Cervical Vertebra Tools               ")
    print("=======================================================")           
        
    ScriptedLoadableModuleWidget.setup(self)
    
    # to access logic class functions and setup global variables
    self.logic = CervicalVertebraToolsLogic()
    # Set default VisSIm location in the user home 
    #TODO: add option user-defined path when installed first time 
    self.vtVars = self.logic.setGlobalVariables(True)
    
    #=================================================================
    #                     Create the GUI interface
    #=================================================================   
    # Create main collapsible Button 
    self.mainCollapsibleBtn = ctk.ctkCollapsibleButton()
    self.mainCollapsibleBtn.setStyleSheet("ctkCollapsibleButton { background-color: DarkSeaGreen  }")
    self.mainCollapsibleBtn.text = "VisSim Cervical Vertebra Tools"
    self.layout.addWidget(self.mainCollapsibleBtn)
    self.mainFormLayout = qt.QFormLayout(self.mainCollapsibleBtn)
  
    # Create input Volume Selector
    self.inputSelectorCoBx = slicer.qMRMLNodeComboBox()
    self.inputSelectorCoBx.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.inputSelectorCoBx.setFixedWidth(200)
    self.inputSelectorCoBx.selectNodeUponCreation = True
    self.inputSelectorCoBx.addEnabled = False
    self.inputSelectorCoBx.removeEnabled = False
    self.inputSelectorCoBx.noneEnabled = False
    self.inputSelectorCoBx.showHidden = False
    self.inputSelectorCoBx.showChildNodeTypes = False
    self.inputSelectorCoBx.setMRMLScene( slicer.mrmlScene )
    self.inputSelectorCoBx.setToolTip("select the input image")
    self.mainFormLayout.addRow("Input image: ", self.inputSelectorCoBx)

    # use qtcombobox
    self.vILbl = qt.QLabel()
    self.vILbl.setText("Which Vertebra? 1-7")        
    self.vILbl.setFixedHeight(20)
    self.vILbl.setFixedWidth(150)
    
    #TODO: include head and shoulders
    self.vtIDCoBx = qt.QComboBox()
    self.vtIDCoBx.addItems(["C1","C2","C3","C4","C5","C6","C7"])
    self.vtIDCoBx.setCurrentIndex(2)
    self.vtIDCoBx.setFixedHeight(20)
    self.vtIDCoBx.setFixedWidth(100)        
    #self.vtIDCoBx.setReadOnly(False) # The point can only be edited by placing a new Fiducial
    # if changed , the default value will change  
    self.vtIDCoBx.connect("currentIndexChanged(int)", self.onVtIDCoBxChange)                  
    #self.mainFormLayout.addRow( self.vILbl,  self.vtIDCoBx )        
      
    # Create a textbox for vertebra location
    # TODO activate input IJK values as well
    Pt = [0,0,0]
    self.inputPointEdt = qt.QLineEdit()
    self.inputPointEdt.setFixedHeight(20)
    self.inputPointEdt.setFixedWidth(100)
    self.inputPointEdt.setText(str(Pt))
    self.inputPointEdt.connect("textChanged(QString)", self.onInputPointEdtChanged)                                  
    #self.inputPointEdt.connect("textEdited(str)", self.onInputPointEdtEdited)                                  

    self.mainFormLayout.addRow( self.inputPointEdt, self.vtIDCoBx)    


    # Add check box for extracting ligaments points 
    self.ligPtsChkBx = qt.QCheckBox("Ligaments points")
    self.ligPtsChkBx.checked = True
    self.ligPtsChkBx.stateChanged.connect(self.onLigPtsChkBxChange)
    self.mainFormLayout.addRow(self.ligPtsChkBx)

    # Create a time label
    self.timeLbl = qt.QLabel("  Time: 00:00")
    self.timeLbl.setFixedWidth(500)   
    self.tmLbl = self.timeLbl
    
    # Create a button to run segmentation
    self.applyBtn = qt.QPushButton("Run")
    self.applyBtn.setFixedHeight(50)
    self.applyBtn.setFixedWidth (150)
    self.applyBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
    self.applyBtn.toolTip = ('How to use:' ' Load an images into Slicer. Pick vertebra locations using the buttons and the Slicer Fiducial tool ')
    self.applyBtn.connect('clicked(bool)', self.onApplyBtnClick)
    self.mainFormLayout.addRow(self.applyBtn, self.timeLbl)
    self.runBtn = self.applyBtn

    # Create a button to display result folder
    self.openResultFolderBtn = qt.QPushButton("open output folder")
    self.openResultFolderBtn.setFixedHeight(20)
    self.openResultFolderBtn.setFixedWidth (150)
    self.openResultFolderBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
    self.openResultFolderBtn.toolTip = ('Scale model and points to 1 mm and translate them to the model center of mass ')
    self.openResultFolderBtn.connect('clicked(bool)', self.onOpenResultFolderBtnClick)
    
    self.mainFormLayout.addRow(self.openResultFolderBtn)

    self.layout.addStretch(1) # Collapsible button is held in place when collapsing/expanding.
    lm = slicer.app.layoutManager();    lm.setLayout(2)

  def cleanup(self):#nothing to do
    pass
  #enddef

  #------------------------------------------------------------------------
  #                        Vertebra Selection
  #------------------------------------------------------------------------
  def onVtIDCoBxChange(self):
      self.inputVolumeNode=self.inputSelectorCoBx.currentNode()
      self.vtID = self.vtIDCoBx.currentIndex + 1   
      self.logic.inputVolumeNode  =  self.inputVolumeNode 
      self.logic.inputFiducialNode = None
      nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
      for f in nodes:
          if ((f.GetName() == self.inputVolumeNode.GetName()+"_vtLocations") ):
             #replace  current 
             print("inputFiducialNode exist")
             self.logic.inputFiducialNode = f  
             newNode= False
            #endif
      #endfor    

      self.logic.setVtID(self.vtID, self.logic.inputVolumeNode , self.logic.inputFiducialNode)
      self.logic.locateVertebra(self.inputVolumeNode, self.vtID, self.inputPointEdt)    
      #self.onInputFiducialBtnClick()
  #enddef

  def onInputPointEdtChanged(self,point):
      self.vtID = self.vtIDCoBx.currentIndex + 1   
      self.logic.setVtIDfromEdt(point,self.vtID)
  #enddef
    
  # extract ligaments points 
  def onLigPtsChkBxChange(self):      
      nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
      self.logic.setLigChk(self.ligPtsChkBx.checked, nodes)
  #enddef
  
  
  def onApplyBtnClick(self):
    self.runBtn.setText("...please wait")
    self.runBtn.setStyleSheet("QPushButton{ background-color: red  }")
    slicer.app.processEvents  
    self.stm=time.time()
    print("time:" + str(self.stm))
    self.timeLbl.setText("                 Time: 00:00")

    self.vtID = self.vtIDCoBx.currentIndex + 1
    self.inputNode = self.inputSelectorCoBx.currentNode()
    pointSelected = self.inputPointEdt.text =="[0,0,0]"
    try: 
       if  (not self.inputNode is None) and (not pointSelected) and (not self.logic.inputFiducialNode is None):
            # create an option to use IJK point or fidicual node
            # inputImage, FiducialPoint, vertebraID, isExteranl ,C7Editbox
            self.logic.run(self.inputNode ,self.logic.inputFiducialNode, self.vtID ,False)

       else:
            print("error in input")
    except Exception as e:
            print("STOPPED: error in input")
            print(e)
    #endtry        
            
    self.etm=time.time()
    tm=self.etm - self.stm
    self.timeLbl.setText("Time: "+str(tm)+"  seconds")
    self.runBtn.setText("Run")
    self.runBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
    slicer.app.processEvents()
  #enddef
       
  def onOpenResultFolderBtnClick(self):
      output = expanduser("~/VisSimTools")   + "/outputs"
      self.logic.openResultsFolder(output)
  #enddef      

#===================================================================
#                           Logic
#===================================================================
class CervicalVertebraToolsLogic(ScriptedLoadableModuleLogic):


  ElastixLogic = Elastix.ElastixLogic()
  ElastixBinFolder = ElastixLogic.getElastixBinDir() 
  vtVars = {}

  # string to boolean converter
  def s2b(self,s):
        return s.lower() in ("yes", "true", "t", "1")
  #enddef
  
  #set global paths and parameters
  def setGlobalVariables(self, isExternalCall):
    print("SpineToolsLogic: initializing global variables:")  
    #define a dictonary
    self.vtVars['vissimPath']           = os.path.join(expanduser("~"),"VisSimTools")
    self.vtVars['elastixBinPath']       = os.path.join(self.ElastixBinFolder, "elastix")
    self.vtVars['transformixBinPath']   =  os.path.join(self.ElastixBinFolder, "transformix")  
    self.vtVars['othersWebLink']        = "https://cloud.uni-koblenz-landau.de/s/yfwcdymS9QfqKc9/download"
    self.vtVars['noOutput']             = " >> /dev/null"
    self.vtVars['outputPath']           = os.path.join(self.vtVars['vissimPath'],"outputs")
    parsPath                            = self.vtVars['vissimPath']  + ",pars,parSpiSeg.txt" 
    self.vtVars['parsPath']             = os.path.join(*parsPath.split(","))
    modelPath                           = self.vtVars['vissimPath']  + ",models,modelCervicalSpine" 
    self.vtVars['modelPath']            = os.path.join(*modelPath.split(","))
    self.vtVars['imgType']              = ".nrrd"
    self.vtVars['vtID']                 = "7"
    vtMethodsegT= [",Default"]
    vtMethodsgT = ["S.seg" ]
    self.vtVars['vtMethodID']           = "0"
    self.vtVars['segT']                 = vtMethodsegT[int(self.vtVars['vtMethodID'])] 
    self.vtVars['sgT']                  = vtMethodsgT[int(self.vtVars['vtMethodID'])] 
    self.vtVars['Styp']                 = "Ct"  
    self.vtVars['vtPtsLigDir']          =  ",PtsLig"
    self.vtVars['vtPtsLigSuff']         =  "Lp"
    modelCropImgLigPtsTmpPath           =   self.vtVars['modelPath'] +"," +self.vtVars['vtPtsLigDir']+","+self.vtVars['Styp']
    self.vtVars['modelCropImgLigPtsTmpPath'] =   os.path.join(*modelCropImgLigPtsTmpPath.split(","))  
    subVarsTemplateFnm                       = self.vtVars['modelPath'] +","+self.vtVars['vtPtsLigDir']+",simPackSubVars.txt"  
    self.vtVars['subVarsTemplateFnm']        =   os.path.join(*subVarsTemplateFnm.split(","))  
    self.vtVars['dispViewTxt']               = "Yellow"
    self.vtVars['dispViewID']                = "7"
    self.vtVars['downSz']                    = "160"
    self.vtVars['winOS']                     = "False"
    self.vtVars['ligChk']                    = "True"
    self.vtVars['hrChk']                     = "True"
    self.vtVars['segNodeCoM']                = "[ 0 , 0 , 0 ]"
    self.vtVars['croppingLength']            = "[ 80 , 80 , 50 ]"   
    self.vtVars['RSxyz']                     = "[ 0.5, 0.5 , 0.5 ]"
    # change the model type from vtk to stl 
    msn=slicer.vtkMRMLModelStorageNode()
    msn.SetDefaultWriteFileExtension('stl')
    slicer.mrmlScene.AddDefaultNode(msn)

    if sys.platform == 'win32':
           self.vtVars['elastixBinPath']       = os.path.join(self.ElastixBinFolder, "elastix.exe")
           self.vtVars['transformixBinPath']   = os.path.join(self.ElastixBinFolder, "transformix.exe")
    #endif   

    self.checkVisSimTools(self.vtVars)
    return self.vtVars
    
  # Check if image is valid
  def hasImageData(self,inputVolumeNode):
    #check input image 
    if not inputVolumeNode:
      logging.debug('hasImageData failed: no input volume node')
      return False
    if inputVolumeNode.GetImageData() is None:
      logging.debug('hasImageData failed: no image data in input volume node')
      return False
    return True
  #enddef

  # set resampling option
  def setLigChk(self,ligChk,nodes):
        self.vtVars['ligChk'] = str(ligChk)
        #self.vtPtsLigTxt     = "Ligaments points"     
        #self.vtVars['vtPtsLigDir']     = "PtsLig" 
        #self.vtVars['vtPtsLigSuff']    = "Lp"
        ligName = "_LigPts_C" # + self.vtVars['vtID'] 
        # Set unvisible 
        for f in nodes:
            if (ligName in f.GetName() ):
               f.GetDisplayNode().SetVisibility(self.s2b(self.vtVars['ligChk']))
               print("Ligaments points detection:" + str(ligChk))
               #break
            #endif
        #endfor
  #enddef


  # set resampling option
  def setHrChk(self,hrChk):
      self.vtVars['hrChk'] = str(hrChk)
  #enddef
  
  def setmrChk(self,mrChk):
      self.mrChk = mrChk
      if mrChk:       
            self.vtVars['Styp']="Mr"
      else:
            self.vtVars['Styp']="Ct"    
      #endif
  #enddef    

  def setVtID(self,idx,inputVolumeNode , inputFiducialNode):
      self.vtVars['vtID']=str(idx)
      print(self.vtVars['vtID']+ " is selected")
      self.inputVolumeNode = inputVolumeNode
      self.inputFiducialNode = inputFiducialNode 
      #remove old points 
        
      # Check if a markup node exists
      newNode = True
      print(self.inputVolumeNode.GetName()+"_vtLocations")
      self.inputFiducialNode = None
      nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
      for f in nodes:
          if ((f.GetName() == self.inputVolumeNode.GetName()+"_vtLocations") ):
             #replace  current 
             print("inputFiducialNode exist")
             self.inputFiducialNode = f  
             newNode= False
            #endif
      #endfor      
      if not (self.inputFiducialNode is None):
         ls = slicer.modules.markups.logic()
         ls.SetActiveListID(self.inputFiducialNode)
         print(ls.GetActiveListID())
         noPts = self.inputFiducialNode.GetNumberOfFiducials() 
         newFid= True
         for j in range (0, noPts):
             if self.inputFiducialNode.GetNthFiducialLabel(j)==("C"+self.vtVars['vtID']) :
                newFid= False 
                print("C"+self.vtVars['vtID'] +" exist, removing old point at: " +str(j))
                #get the new location
                self.inputFiducialNode.RemoveMarkup(j)      
             #endif
         #endfor
      else:
         print("inputFiducialNode does not exist")

         
      #endif  
      
  #enddef
  
  # check if vertebra location is available
  def setVtIDfromEdt(self,point,vtID):
        # no external call
        #print("no external call, point= " + point)
        self.vtVars['vtID'] = str(vtID)
        isExternalCall = False
        print("point changed,  " + str(vtID) + " is selected")      
        #TODO: add option to use point from text for cropping   
        return isExternalCall     
  #enddef

  def locateVertebra(self, inputVolumeNode, vtID,  inputPointEdt):
        # redefine to be used in the logic class.
        self.inputVolumeNode   = inputVolumeNode
        self.inputPointEdt     = inputPointEdt  
        self.vtVars['vtID']= str(vtID)
        print(" ..... getting vertebra location in the input image")  
        # Reset global point label
        self.inputPoint = [0,0,0]
        self.inputPointEdt.setText(str(self.inputPoint))
        # Check if a volume is selected
        if not self.inputVolumeNode:
           print >> sys.stderr, "You need to pick a input volume first before locating vertebra."
           return -1
        #endif
        
        #  Display Sagittal during locating the vertebra
        disp_logic = slicer.app.layoutManager().sliceWidget(self.vtVars['dispViewTxt']).sliceLogic()
        disp_cn = disp_logic.GetSliceCompositeNode()
        disp_cn.SetBackgroundVolumeID(self.inputVolumeNode.GetID())
        lm = slicer.app.layoutManager()
        #slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpGreenSliceView 8
        #slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpYellowSliceView 7       
        lm.setLayout(int(self.vtVars['dispViewID']))
        # Fit slice to window
        sliceNodes = slicer.util.getNodes('vtkMRMLSliceNode*')
        layoutManager = slicer.app.layoutManager()
        for sliceNode in sliceNodes.values():
            sliceWidget = layoutManager.sliceWidget(sliceNode.GetLayoutName())
            if sliceWidget:
                sliceWidget.sliceLogic().FitSliceToAll()
            #endif
        #endfor
       
        #TODO: remove this part, it is already implemented in setVtID
        # Check if a markup node exists
        newNode = True
        nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
        for f in nodes:
            if ((f.GetName() == inputVolumeNode.GetName()+"_vtLocations") ):
                #replace  current 
                self.inputFiducialNode = f  
                newNode= False
            #endif
        #endfor
        
        #this is implemented in setVtID
        if newNode:
           self.inputFiducialNode = slicer.vtkMRMLMarkupsFiducialNode()
           self.inputFiducialNode.SetName(inputVolumeNode.GetName()+"_vtLocations")
           slicer.mrmlScene.AddNode(self.inputFiducialNode)
        #endif   
        self.inputFiducialNode.GetDisplayNode().SetTextScale(2)
        self.inputFiducialNode.GetDisplayNode().SetSelectedColor(1,0,0)           

        # Start Fiducial Placement Mode in Slicer for one node only
        placeModePersistance = 0 # one node only
        slicer.modules.markups.logic().StartPlaceMode(placeModePersistance)

        # Observe scene for updates
        self.addObs = self.inputFiducialNode.AddObserver(self.inputFiducialNode.MarkupAddedEvent,   self.onInputFiducialNodeMarkupAddedEvent)
        self.modObs = self.inputFiducialNode.AddObserver(self.inputFiducialNode.PointModifiedEvent, self.onInputFiducialNodePointModifiedEvent)
        self.rmvObs = self.inputFiducialNode.AddObserver(self.inputFiducialNode.MarkupRemovedEvent, self.onInputFiducialNodeMarkupRemovedEvent)

  #enddef

  def onInputFiducialNodeMarkupAddedEvent(self, caller, event):
        # it seems this action happened after adding new fiducial
        print("Fiducial adding event!")
        #remove previous observer
        caller.RemoveObserver(self.addObs)
        noPts = caller.GetNumberOfFiducials() 
        rasPt = [0,0,0] 
        caller.SetNthFiducialLabel(noPts-1, "C"+self.vtVars['vtID'])
        caller.GetNthFiducialPosition(noPts-1,rasPt)
        self.inputPoint = self.ptRAS2IJK(caller, noPts-1, self.inputVolumeNode)
        self.inputPointEdt.setText(str(self.inputPoint))
        #self.inputPointEdt.setText(str(rasPt))
        print(" ..... vertebra location RAS: " + str(rasPt))  
        print(" ..... vertebra location in the fixed image set to: " + str(self.inputPoint))
  #enddef
  
  
  #--------------------------------------------------------------------------------------------
  #    RAS to  IJK Event
  #--------------------------------------------------------------------------------------------
  # The fiducial point saved in RAS, we need to convert to IJK
  #  more info in our wiki 
  def onInputFiducialNodePointModifiedEvent(self, caller, event):
      print("Fiducial modified event!")
      #caller.RemoveObserver(self.modObs)
      # get the new IJK position and display it
      rasPt = [0,0,0] 
      i = caller.GetAttribute('Markups.MovingMarkupIndex')
      if not (i is None):
         i=int(i)
         caller.GetNthFiducialPosition(i,rasPt)
         #self.inputPointEdt.setText(str(rasPt))
         self.inputPoint = self.ptRAS2IJK(caller, i, self.inputVolumeNode)
         self.inputPointEdt.setText(str(self.inputPoint))

      #endif   
  #enddef

  def onInputFiducialNodeMarkupRemovedEvent(self, caller, event):
      print("Fiducial removed event!")
      caller.RemoveObserver(self.rmvObs)
      #i = caller.GetNumberOfFiducials()-1
      #print("number of rmaining fiducials: " + str(i))
      
  def openResultsFolder(self,outputPath):
      if outputPath is None:
         outputPath =self.vtVars['outputPath'] 
      #endif
      if sys.platform == 'darwin':
        os.system('open ' + outputPath)
      elif sys.platform == 'linux2':
        os.system('xdg-open ' + outputPath)
      elif sys.platform == 'win32':
        s = 'explorer  ' + outputPath
        s = s.replace('/', '\\')        
        os.system(s)
      #endif
  #enddef
  
  #--------------------------------------------------------------------------------------------
  #                        run elastix
  #--------------------------------------------------------------------------------------------      
  def runElastix(self, elastixBinPath, fixed, moving, output, parameters, verbose, line):
       print ("************  Compute the Transform **********************")
       Cmd =elastixBinPath + " -f " + fixed +" -m "+ moving  + " -out " + output  + " -p " + parameters + verbose
       print("Executing: " + Cmd)
       #cTI=os.system(cmd)
       si = None
       if sys.platform == 'win32':
         si = subprocess.STARTUPINFO()
         si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
       #endif  
       cTI = subprocess.call(Cmd , shell = (sys.platform == 'linux2') , startupinfo=si )
       errStr="elastix error at line"+ line +", check the log files"
       self.chkElxER(cTI,errStr) # Check if errors happen during elastix execution
       return cTI
  #enddef

  #--------------------------------------------------------------------------------------------
  #                        run transformix
  #--------------------------------------------------------------------------------------------      
  def runTransformix(self,transformixBinPath, img, output, parameters, verbose, line):
      print ("************  Apply transform **********************")
      si = None
      if sys.platform == 'win32':
         si = subprocess.STARTUPINFO()
         si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
      #endif  		
      # Apply the transformation to the segmentation image:
      Cmd = transformixBinPath + " -in " + img + " -out " + output  + " -tp " + parameters +" -def all " +verbose
      print("Executing... " + str(Cmd))
      #cTS=os.system(Cmd)
      cTS = subprocess.call(Cmd , shell = (sys.platform == 'linux2') , startupinfo=si )
      errStr="transformix error at line"+ line +", check the log files"
      self.chkElxER(cTS,errStr) # Check if errors happen during elastix execution
      return cTS
  #enddef

  #--------------------------------------------------------------------------------------------
  #                       Cropping Process  
  #--------------------------------------------------------------------------------------------
  # Using the location as a center point, we cropp around it using the defined cropLength 
  def runCropping(self, inputVolume, point, vtIDt,croppingLengthT, samplingLengthT, hrChkT):
        print("================= Begin cropping ... =====================")
        # Create a temporary node as workaround for bad path or filename 
        #TODO: create a temp folder and remove temp node before display
        vtID = int(vtIDt)
        croppingLength =   self.t2v(croppingLengthT)
        samplingLength =   self.t2v(samplingLengthT)
        hrChk          =   self.s2b(hrChkT)
        
        print("Vertebra "+str(vtID)+ " location: " + str(point) + "   cropping length: " + str(croppingLength) )
        #Remove old cropping node
        nodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
        for f in nodes:
            if ("inputimage_crop" in f.GetName()):
                slicer.mrmlScene.RemoveNode(f )
            #endif
        #endfor
        # resampling spacing  
        self.RSx= samplingLength[0] ; self.RSy=samplingLength[1];     self.RSz= samplingLength[2]

        #Get input image information 
        spacing = inputVolume.GetSpacing()
        imgData = inputVolume.GetImageData()
        dimensions = imgData.GetDimensions()
        
        # compute cropping bounds from image information and cropping parameters
        croppingBounds = [[0,0,0],[0,0,0]];   size = [0,0,0];    lower = [0,0,0] ;     upper = [0,0,0]
        for i in range(0,3):
            size[i] = int((croppingLength[i]/spacing[i])/2)
            lower[i] = int(point[i]) - int(size[i])
            upper[i] = dimensions[i] - int(point[i]+size[i])
            # Check if calculated boundaries exceed image dimensions
            if lower[i] < 0:
                    lower[i] = 0
            #endif        
            if upper[i] > dimensions[i]:
                   upper[i] = int(dimensions[i])
            #endif
        #endfor   
        croppingBounds = [lower,upper]
        # Call SimpleITK CropImageFilter
        print("Cropping with " + str(croppingBounds[0]) + " and " + str(croppingBounds[1]) + ".")
        #self.inputCropPath = "bla blaa blaaa"

        inputImage = sitkUtils.PullVolumeFromSlicer(inputVolume.GetID())
        cropper = sitkUtils.sitk.CropImageFilter()
        croppedImage = cropper.Execute(inputImage, croppingBounds[0], croppingBounds[1])          
        nodeName = str(inputVolume.GetName()) +"_C"+str(vtID) +"_crop"
        inputCropPath = os.path.splitext(inputVolume.GetStorageNode().GetFileName())[0] +"_C"+str(vtID) +"_crop.nrrd"
        # Make a node with cropped image 
        sitkUtils.PushVolumeToSlicer(croppedImage, None, nodeName , 'vtkMRMLScalarVolumeNode' )
        crNodes = slicer.util.getNodesByClass("vtkMRMLScalarVolumeNode")
        for f in crNodes:
            if nodeName in f.GetName():
                 f.SetName(nodeName) 
                 break         
            #endif
        #endfor  
        croppedNode = slicer.util.getNode(nodeName)
        inputCropPath = self.vtVars['vissimPath']+",inputImage"  +"_C"+str(vtID) +"_crop.nrrd"                                  
        print(inputCropPath.split(","))
        inputCropPath = os.path.join(*inputCropPath.split(","))                         
        print("cropped:     "+inputCropPath)
        slicer.util.saveNode( croppedNode, inputCropPath)

        #-------------------------------------------------------
        # Resampling: this produces better looking models  
        #-------------------------------------------------------
        if hrChk:
           #Run slicer cli module: resample scalar volume
           #inputCropIsoPath = os.path.splitext(inputVolume.GetStorageNode().GetFileName())[0] +"_C"+str(vtID) +"_crop_iso.nrrd"  
           inputCropIsoPath = self.vtVars['vissimPath']+",inputImage"  +"_C"+str(vtID) +"_crop_iso.nrrd"                                  
           inputCropIsoPath = os.path.join(*inputCropIsoPath.split(","))                         
           print("iso cropped: "+inputCropIsoPath)
           resampleSpacing = " ["+ str(self.RSx) + "," + str(self.RSy) + "," + str(self.RSz) + "] "
           #TODO: Get Slicer PATH
           SlicerPath =os.path.abspath(os.path.join(os.path.abspath(os.path.join(os.sys.executable, os.pardir)), os.pardir))
           SlicerBinPath =   os.path.join(SlicerPath,"Slicer")  
           ResampleBinPath = SlicerPath +",lib,Slicer-4.10,cli-modules,ResampleScalarVolume"
           ResampleBinPath =  os.path.join(*ResampleBinPath.split(",")) 		   
           resamplingCommand = SlicerBinPath + " --launch " + ResampleBinPath   
           si = None 
           if sys.platform == 'win32':
              #note: in windows, no need to use --launch
              resamplingCommand = ResampleBinPath + ".exe"
              print(os.path.getsize(resamplingCommand))
              si = subprocess.STARTUPINFO()
              si.dwFlags |= subprocess.STARTF_USESHOWWINDOW
		   #endif

           cmdPars = " -i linear -s "+ resampleSpacing + inputCropPath +" "+inputCropIsoPath  
           Cmd = resamplingCommand  + cmdPars
           print("Executing ... "+Cmd)
           #os.system(cmd)
           subprocess.call(Cmd , shell = (sys.platform == 'linux2') , startupinfo=si )
           #inputCropPath = inputCropIsoPath
           print(" Cropping and resampling are done !!! ")
        #endif

        #inputCropPath    = inputCropPath.strip("'")
        print(" Cropped image is saved in : [%s]" % inputCropPath)
        print(" Cropping is done !!! ")          

        return inputCropPath
  
  #===========================================================================================
  #                       Segmentation Process 
  #--------------------------------------------------------------------------------------------        
  # This method perform the atlas segementation steps
  def run(self, inputVolumeNode, inputFiducialNode, vtID, isExternalCall):
      #to be used fromoutside we need to do:
      # import CervicalVertebraTools
        """
        Run the actual algorithm
        """
        logging.info('Processing started')

        if isExternalCall:
           print(" External call!")
           # this functionshould be changed in case of external call
           self.vtVars = self.setGlobalVariables(isExternalCall)
        else:
           print(" No external call!")
        
        vtID = int(vtID)
        print (vtID)    
   
        # set the correct models path:
		
        self.modelCropImgLigPtsPath = self.vtVars['modelCropImgLigPtsTmpPath'] +str(vtID)+self.vtVars['vtPtsLigSuff']+".fcsv"

        # we need to run this again in case of external call
        #endif
        modelCropPath       =   self.vtVars['modelPath']+  ',Default'+",Mdl" + self.vtVars['Styp']+ str(vtID)        +self.vtVars['imgType'] 
        modelCropPath       =   os.path.join(*modelCropPath.split(","))
        print(modelCropPath)
        modelCropSegPath    = self.vtVars['modelPath'] + self.vtVars['segT']+",Mdl"+self.vtVars['Styp']+ str(vtID)+ self.vtVars['sgT'] +self.vtVars['imgType'] 
        modelCropSegPath    =   os.path.join(*modelCropSegPath.split(","))

        self.hasImageData(inputVolumeNode)

        if not os.path.exists(self.vtVars['outputPath']):
           os.mkdir(self.vtVars['outputPath'])      
        else:
           #only for this vertebra
           self.removeOtputsFolderContents()
        #endif
        # results paths        
        resTransPath  = os.path.join(self.vtVars['outputPath'] ,"TransformParameters.0.txt")
        resOldDefPath = os.path.join(self.vtVars['outputPath'] , "deformationField"+self.vtVars['imgType'])
        resDefPath    = os.path.join(self.vtVars['outputPath'] , inputVolumeNode.GetName()+"C"+str(vtID)+"_dFld"+self.vtVars['imgType'])
        
        #remove old result files:
        if os.path.isfile(resOldDefPath):
           os.remove(resOldDefPath) 
        if os.path.isfile(resDefPath):
           os.remove(resDefPath)
        
        # check if the model is found
        if not isfile(modelCropPath): 
            print >> sys.stderr, "ERROR: model is not found"            
            print("modelPath: " + modelCropPath)
            return False
        # endif

        inputPath = inputVolumeNode.GetStorageNode().GetFileName()
        inputFnm  = basename(os.path.splitext(inputPath)[0])    
        
        #remove old results
        resultSegNodeName    = inputVolumeNode.GetName()    + "_Seg_C"      +str(vtID)           
        resultLigPtsNodeName = inputVolumeNode.GetName()    + "_LigPts_C"   +str(vtID)
        resultTransformNodeName =  inputVolumeNode.GetName()+ "_Transform_C"+str(vtID)
            
        # Get IJK point from the fiducial to use in cropping  
        newFid= True
        for j in range (inputFiducialNode.GetNumberOfFiducials() ):
             if inputFiducialNode.GetNthFiducialLabel(j)==("C"+str(vtID)) :
                break
              #endif
        #endfor
        inputPoint = self.ptRAS2IJK(inputFiducialNode,j,inputVolumeNode)

        #Remove old results
        for node in slicer.util.getNodes():
            if ( resultLigPtsNodeName    == node): slicer.mrmlScene.RemoveNode(node) #endif
            if ( resultSegNodeName       == node): slicer.mrmlScene.RemoveNode(node) #endif
            if ( resultTransformNodeName == node): slicer.mrmlScene.RemoveNode(node) #endif
        #endfor    
         
        # TODO: add better condition
        if  np.sum(inputPoint)== 0 :
            print("Error: select vertebra point")
            return False
        #endif  

        print ("************  Cropping  **********************")
        self.vtVars['intputCropPath'] = self.runCropping(inputVolumeNode, inputPoint, vtID ,self.vtVars['croppingLength'], self.vtVars['RSxyz'],  self.vtVars['hrChk'] )                     
        [success, croppedNode] = slicer.util.loadVolume(self.vtVars['intputCropPath'], returnNode=True)
        print(self.vtVars['intputCropPath'])
        croppedNode.SetName(inputVolumeNode.GetName()+"_C"+str(vtID)+"Crop")                                                        
        print ("************  Register model to cropped input image **********************")
        cTI = self.runElastix(self.vtVars['elastixBinPath'],self.vtVars['intputCropPath'],  modelCropPath, self.vtVars['outputPath'], self.vtVars['parsPath'], self.vtVars['noOutput'], "554")
       
		#genrates deformation field 
        cTR = self.runTransformix(self.vtVars['transformixBinPath'],modelCropPath, self.vtVars['outputPath'], resTransPath, self.vtVars['noOutput'], "556")
        # rename fthe file:
        os.rename(resOldDefPath,resDefPath)
        
        print ("************  Load deformation field Transform  **********************")
        [success, vtTransformNode] = slicer.util.loadTransform(resDefPath, returnNode = True)
        vtTransformNode.SetName(resultTransformNodeName)

        print ("************  Transform segmentation  **********************")
        print(modelCropSegPath)
        [success, vtResultSegNode] = slicer.util.loadSegmentation(modelCropSegPath, returnNode = True)
        vtResultSegNode.SetName(resultSegNodeName)
        vtResultSegNode.SetAndObserveTransformNodeID(vtTransformNode.GetID()) # movingAllMarkupNode should be loaded, the file contains all points
        slicer.vtkSlicerTransformLogic().hardenTransform(vtResultSegNode) # apply the transform
        vtResultSegNode.CreateClosedSurfaceRepresentation() 
                     
        if self.s2b(self.vtVars['ligChk']):            
           print ("************  Transform Ligaments Points **********************")
           modelCropImgLigPtsPath = self.vtVars['modelPath'] +","+self.vtVars['vtPtsLigDir']+","+self.vtVars['Styp']+ str(vtID)+self.vtVars['vtPtsLigSuff']+".fcsv"
           modelCropImgLigPtsPath = os.path.join(*modelCropImgLigPtsPath.split(","))
           print(self.modelCropImgLigPtsPath)
           [success, vtResultLigPtsNode] = slicer.util.loadMarkupsFiducialList  (modelCropImgLigPtsPath, returnNode = True)
           vtResultLigPtsNode.GetDisplayNode().SetTextScale(1)
           vtResultLigPtsNode.GetDisplayNode().SetSelectedColor(1,0,0)           
           vtResultLigPtsNode.SetName(resultLigPtsNodeName)

           vtResultLigPtsNode.SetAndObserveTransformNodeID(vtTransformNode.GetID()) # movingAllMarkupNode should be loaded, the file contains all points
           slicer.vtkSlicerTransformLogic().hardenTransform(vtResultLigPtsNode) # apply the transform
           # needed in extract scaled model
           self.vtResultLigPtsNode = vtResultLigPtsNode       
        #endif 
        # Display the result if no error
        # Clear vertebra location labels
        if  (cTI==0) and (cTR==0):
             # change the model type from vtk to stl 
             msn=slicer.vtkMRMLModelStorageNode()
             msn.SetDefaultWriteFileExtension('stl')
             slicer.mrmlScene.AddDefaultNode(msn)

             print("display result")
             self.dispVolume(inputVolumeNode)
        
             print("get vertebra information")
             tableName =  inputVolumeNode.GetName()+"_tbl"
             # create only if it does not exist
             try:
                 resultsTableNode =  slicer.util.getNode(tableName)
             except:             
                 resultsTableNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode")
                 resultsTableNode.SetName(tableName)
                 resultsTableNode.AddEmptyRow()    
                 resultsTableNode.GetTable().GetColumn(0).SetName("Vertebra")
                 resultsTableNode.AddColumn()
                 resultsTableNode.GetTable().GetColumn(1).SetName("Volume mm3")
                 resultsTableNode.AddColumn()
                 resultsTableNode.GetTable().GetColumn(2).SetName("CoM X")
                 resultsTableNode.AddColumn()
                 resultsTableNode.GetTable().GetColumn(3).SetName("CoM Y")
                 resultsTableNode.AddColumn()
                 resultsTableNode.GetTable().GetColumn(4).SetName("CoM Z")
             #endif
                     
             resultsTableNode = self. getVertebraInfo( vtResultSegNode, croppedNode, vtID, resultsTableNode)
             resultsTableNode.RemoveRow(resultsTableNode.GetNumberOfRows())   
             
             #slicer.app.layoutManager().setLayout( slicer.modules.tables.logic().GetLayoutWithTable(slicer.app.layoutManager().layout))
             slicer.app.applicationLogic().GetSelectionNode().SetActiveTableID(resultsTableNode.GetID())
             slicer.app.applicationLogic().PropagateTableSelection()
             slicer.app.layoutManager().setLayout(1)            
        else:
            print("error happened during segmentation ")
        #endif
        
        #Remove temporary nodes:
        os.remove(self.vtVars['intputCropPath'])     
        if not  isExternalCall:     
            slicer.mrmlScene.RemoveNode(croppedNode )
        #endif
        # needed in extract scaled model
        self.vtResultSegNode = vtResultSegNode
        
        self.removeOldFiles(self.vtVars['outputPath'])                            
        print("================= vertebra analysis is complete  =====================")
        logging.info('Processing completed')
        return True
    #enddef
 
  #--------------------------------------------------------------------------------------------
  #                       Check Elastix error
  #--------------------------------------------------------------------------------------------
  # This method checks if errors happen during elastix execution
  def chkElxER(self,c, s):
        if c>0:
           #qt.QMessageBox.critical(slicer.util.mainWindow(),'segmentation', s)
           print(s)  
           return False
        else: 
            print("done !!!")
        #endif
 #enddef 
  
  # segmenteditor effect on the resulted segmentations 
  # this function is called by functions like doSmoothing and doMargining
  def getSegmentationEditor(self,segNode,masterNode):
            # Create segment editor to get access to effects
           segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
           segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
           segmentEditorNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentEditorNode")
           segmentEditorWidget.setMRMLSegmentEditorNode(segmentEditorNode)
           segmentEditorWidget.setSegmentationNode(segNode)
           segmentEditorWidget.setMasterVolumeNode(masterNode)
           return segmentEditorWidget, segmentEditorNode

  # smoothing segmentation effect    
  #TODO: add more options   
  def runSmoothing(self,segNode,masterNode,KernelSizeMm):       
           # Smoothing
           [segEditorW,segEditorN]= self.getSegmentationEditor(segNode,masterNode)               
           for i in range (0,segNode.GetSegmentation().GetNumberOfSegments()):
               segmentID = segNode.GetSegmentation().GetNthSegmentID(i)
               segEditorW.setActiveEffectByName("Smoothing")
               segEditorW.setCurrentSegmentID(segmentID)
               effect = segEditorW.activeEffect()
               effect.setParameter("SmoothingMethod", "MEDIAN")
               effect.setParameter("KernelSizeMm", KernelSizeMm)
               effect.self().onApply()
           #endfor
           # Clean up
           segEditorW = None
           slicer.mrmlScene.RemoveNode(segEditorN)
  #enddef

  # Margin segmentation effect
  # MarginSizeMm>0 Grow, else Shrink          
  def runMargining(self,segNode,masterNode,MarginSizeMm):
           #Dilation and Eroding
           [segEditorW,segEditorN]= self.getSegmentationEditor(segNode,masterNode)               
           for i in range (0,segNode.GetSegmentation().GetNumberOfSegments()):
               segmentID = segNode.GetSegmentation().GetNthSegmentID(i)
               segEditorW.setActiveEffectByName("Margin")
               segEditorW.setCurrentSegmentID(segmentID)
               effect = segEditorW.activeEffect()
               effect.setParameter("MarginSizeMm", MarginSizeMm) 
               effect.self().onApply()
           #endfor
           # Clean up
           segEditorW = None
           slicer.mrmlScene.RemoveNode(segEditorN)
  #enddef

  #--------------------------------------------------------------------------------------------
  #                        Calculate center of mass and volume of a vertebra
  #--------------------------------------------------------------------------------------------
  def getVertebraInfo(self, segNode, masterNode, vtID, resultsTableNode):
        segStatLogic = SegmentStatistics.SegmentStatisticsLogic()
        segStatLogic.getParameterNode().SetParameter("Segmentation", segNode.GetID())
        segStatLogic.getParameterNode().SetParameter("ScalarVolume", masterNode.GetID())
        segStatLogic.getParameterNode().SetParameter("LabelmapSegmentStatisticsPlugin.enabled","False")
        segStatLogic.getParameterNode().SetParameter("ScalarVolumeSegmentStatisticsPlugin.voxel_count.enabled","False")
        segStatLogic.computeStatistics()
        #remove old rows for same vertebra
        for i in range (resultsTableNode.GetNumberOfRows()):
            if "C"+str(vtID) == resultsTableNode.GetCellText(i,0):
               print( "C"+str(vtID) + " tale row exists, old values will be removed.") 
               resultsTableNode.RemoveRow(i)   
        #endif
         
        #export to temporary table then copy to the final table
        resultsTmpTableNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode")
        resultsTmpTableNode.SetName("tmpTable")
        segStatLogic.exportToTable(resultsTmpTableNode)
        idx = resultsTableNode.GetNumberOfRows()-1
        resultsTableNode.AddEmptyRow()    
        resultsTableNode.SetCellText(idx,0,"C"+str(vtID))                          # vertebra
        resultsTableNode.SetCellText(idx,1, resultsTmpTableNode.GetCellText(0,1) ) # volume size
        # if this is C7 compute center of mass        
        if (vtID ==7) and(self.vtVars['vtMethodID']== "0"): # for testing
            segID = segNode.GetSegmentation().GetSegmentIdBySegmentName("C"+str(vtID))
            modelNode = segNode.GetClosedSurfaceRepresentation(segID)
            com = vtk.vtkCenterOfMass(); com.SetInputData(modelNode);   com.Update()
            segNodeCoM = com.GetCenter()
            resultsTableNode.SetCellText(idx,2,str(segNodeCoM[0])) 
            resultsTableNode.SetCellText(idx,3,str(segNodeCoM[1])) 
            resultsTableNode.SetCellText(idx,4,str(segNodeCoM[2]))
            self.vtVars['segNodeCoM']=str(segNodeCoM)
        #endif
        slicer.mrmlScene.RemoveNode(resultsTmpTableNode)
        print("Measurements are completed !!! ")
        return resultsTableNode
  #enddef


  #convert dictonary text to vector
  def t2v(self,txt):
      vector = [0,0,0]
      print(txt)
      t = txt.strip("]").strip("[").strip("(").strip(")")
      t = t.split(",")
      for i in range(3):
          vector[i] =float(t[i])
      return vector 
  #enddef
   
  #enddef
        
  def dispVolume(self,inputVolumeNode):
        lm = slicer.app.layoutManager();    lm.setLayout(4)
        r_logic = lm.sliceWidget("Red").sliceLogic()
        r_cn = r_logic.GetSliceCompositeNode()
        r_cn.SetBackgroundVolumeID(inputVolumeNode.GetID())
        y_logic = lm.sliceWidget("Yellow").sliceLogic()
        y_cn = y_logic.GetSliceCompositeNode()
        y_cn.SetBackgroundVolumeID(inputVolumeNode.GetID())
        g_logic = lm.sliceWidget("Green").sliceLogic()
        g_cn = g_logic.GetSliceCompositeNode()
        g_cn.SetBackgroundVolumeID(inputVolumeNode.GetID())

        #center 3D view and zoom in 3 times
        v3DDWidget = lm.threeDWidget(0)
        v3DDWidgetV = v3DDWidget.threeDView()
        v3DDWidgetV.resetFocalPoint() 
        v3DDWidgetV.zoomFactor =3
        v3DDWidgetV.zoomIn()
        v3DDWidgetV.zoomFactor =0.05 # back to default value

#------------------------------------------------------
#                  IJK to RAS  
#------------------------------------------------------
# This function convert an IJK point to RAS point 
#  input:  a point vector and volume node
#  output: a point vector                     
  def ptIJK2RAS(self,ptIJK,inputImg):
        #TODO: add option for printing                   
        # create a IJK2RAs transformation matrix 
        ijk2rasM = vtk.vtkMatrix4x4()
        inputImg.GetIJKToRASMatrix(ijk2rasM)
        ptRAS=np.zeros((len(ptIJK),3))
        ijk= ptIJK
        # create a 4 elements array to get the converted values
        ijkv=[ijk[0],ijk[1],ijk[2],1]             
        rasPt=ijk2rasM.MultiplyDoublePoint(ijkv)
        ptRAS=rasPt[0:3]
        #print("IJK= " + str(ptIJK)+ "   RAS= " + str(ptRAS))
        return  ptRAS       

#------------------------------------------------------
#                 RAS  to IJK 
#------------------------------------------------------   
# This function convert RAS ro an IJK point 
#  input:  a point vector and volume node
#  output: a point vector                     
  def ptRAS2IJK(self,ptRAS,i,inputImg): 
        #TODO: add option for printing                   
        # create a RAS2IJK transformation matrix 
        ras2ijkM = vtk.vtkMatrix4x4()
        inputImg.GetRASToIJKMatrix(ras2ijkM)       
        ras=[0,0,0]
        #i = int(i)
        print(type(i))
        ptRAS.GetNthFiducialPosition(i,ras)        
        # create a 4 elements array to get the converted values
        rasv=[ras[0],ras[1],ras[2],1]             
        ptIJKf=np.zeros(3);
        ijkPt=ras2ijkM.MultiplyPoint(rasv)
        ptIJKf[0]=ijkPt[0];ptIJKf [1]=ijkPt[1];ptIJKf [2]=ijkPt[2];
        ptIJK = ptIJKf.astype(np.int64)
        #print("RAS= " + str(ras)+ "   IJK= " + str(ptIJK))
        return  ptIJK       
 
  def removeOldFiles(self, outputPath):
      #remove old files if exist
      try:    
          if os.path.isdir(outputPath.strip()): 
             print("removing old output folder!")
             #shutil.rmtree(self.vtVars['outputPath'])
             resfiles = os.listdir(outputPath) 
             oPath= outputPath +"/"
             for fnm in resfiles:
                 print(fnm)
                 if "IterationInfo" in fnm:
                     os.remove(os.path.join(oPath, fnm))
                 elif  "result" in fnm:
                     os.remove(os.path.join(oPath, fnm))
                 elif  ".log" in fnm:
                     os.remove(os.path.join(oPath, fnm))
                 elif  "TransformParameters" in fnm:
                     os.remove(os.path.join(oPath, fnm))
                 #endif
             #endfor                        
           #endif   
      except:
          print("nothing to delete ...")
  #enddef

  def checkVisSimTools(self,vtVars ):
        # TODO: optimise this part to download only the missing files        
        # Check if elastix exist or download it 
        print(" Defaults paths: " + vtVars['vissimPath'])
        print("      VisSimTools folder: " + vtVars['vissimPath'])
        if isfile(vtVars['elastixBinPath'].strip()): 
           print("elastix binaries are found in " + vtVars['elastixBinPath'] )
        else: 
            print("elastix binaries are missing, trying to download ... ")
            self.msgBox("elastix binaries are missing!")
        #endif
        # check if other files exist
        if isfile(vtVars['parsPath'].strip()): 
           print("Other files are found !" )
           print("  Parameter file: " + vtVars['parsPath'])
           print("  Output folder : " + vtVars['outputPath'])            
           #print("  Cropping Length: " + str(croppingLength))           
        else: 
            print("Other files are  missing, trying to download ... ")
            self.msgBox("important files are missing and will be downloaded!")
            try:                               
                print("Downloading VisSimTools others ...")
                vissimZip = expanduser("~/VisSimToolsTmp.zip")
                with open(vissimZip ,'wb') as f:
                     uFile = urllib2.urlopen(vtVars['othersWebLink'])              
                     chunk = 10024096
                     while 1:
                           data = uFile.read(chunk)
                           f.write(data)                   
                           if not data:
                              f.close()                               
                              print "done!"
                              break
                           #endIf
                           print "Reading ...  %s bytes"%len(data) 
                     #endWhile                               
                print("Extracting to user home ")
                zip_ref = zipfile.ZipFile(vissimZip, 'r')
                zip_ref.extractall(expanduser("~/"))
                zip_ref.close()  
                #remove the downloaded zip file     
                os.remove(vissimZip)   
                # change permission of bin folder for Linux
                if int(self.vtVars['downSz'])==0:   
                   print("Making binaries executable for Linux ")
                   md=  stat.S_IRWXU | stat.S_IRGRP | stat.S_IROTH |stat.S_IXGRP |stat.S_IXOTH
                   os.chmod(vtVars['elastixBinPath'].strip()    ,  md)
                   os.chmod(vtVars['transformixBinPath'].strip(),  md)
                   os.chmod(vtVars['elastixBinPath'].strip(),  md)
                #endif 
                msg.setInformativeText("VisSimTools folder is downloaded and ready to use!")
                msg.exec_()                      
                                          
            except Exception as e:
                  print("Error: can not download and extract VisSimTools ...")
                  print(e)   
                  return -1
            #end try-except  
  #enddef
  
  def msgBox(self,txt):
      #msg = qt.QMessageBox()
      #msg.setIcon(qt.QMessageBox.Information)
      #msg.setText("information:")
      #msg.setInformativeText(txt)
      #msg.setWindowTitle("VisSimTools")
      #msg.exec_()
      print(txt)
  #enddef

  
  def removeOtputsFolderContents(self):
      try:
          for file in os.listdir(self.vtVars['outputPath']):
              filePath = os.path.join(self.vtVars['outputPath'], file)
              if os.path.isfile(filePath):
                 os.remove(filePath)
             #endif
          #endfor        			
      except Exception as e:
            print("nothing to delete ...")
            print(e)
       #endtry 
  #enddefr 
              
                                     
#===================================================================
#                           Test
#===================================================================
class CervicalVertebraToolsTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  logic = CervicalVertebraToolsLogic()

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)
    self.vtVars    = self.logic.setGlobalVariables(True)
    self.logic.removeOtputsFolderContents()

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.testSlicerCervicalVertebraTools()

  def testSlicerCervicalVertebraTools(self):

    self.delayDisplay("Starting the test")
    self.stm=time.time()
    print("time:" + str(self.stm))
    
    # to get the links from datastore open http://slicer.kitware.com/midas3/community/23 then select a file and click share to get
    # the download link

    tstImg = 1 # 1 =CT, 2 = MR 
    if tstImg ==1:
       imgWebLink = "https://cloud.uni-koblenz-landau.de/s/Mb6JHLdWw5MEPB2/download"
       fnm = os.path.join(*(self.logic.vtVars['outputPath'] +",imgCT"+self.logic.vtVars['imgType']).split(","))
       c7p = [ 0.390  , -26.445 ,-115.821 ]
    else:
       imgWebLink = "https://cloud.uni-koblenz-landau.de/s/ieyDfHpCjHNpZXi/download"
       fnm = os.path.join(*(self.logic.vtVars['outputPath'] +",imgMR"+self.logic.vtVars['imgType']).split(","))
       c7p = [ -0.408  ,   3.439 ,   54.432 ]
    #endif



    # don't download if already downloaded                       
    if not os.path.exists(fnm):
       try:         
           print("Downloading vertebra sample image ...")
           import urllib
           urllib.urlretrieve (imgWebLink ,fnm )
       except Exception as e:
              print("Error: can not download sample file  ...")
              print(e)   
              return -1
       #end try-except 
    #endif
    # remove old nodes 
    nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
    for f in nodes:
       if ((f.GetName() == "VertebraLocationPoint") ):
           slicer.mrmlScene.RemoveNode(f)
       #endif
    #endfor 
    nodes = slicer.util.getNodesByClass('vtkScalarVolumeNode')
    for f in nodes:
       if (f.GetName() ==  "imgA"):
           slicer.mrmlScene.RemoveNode(f)
       #endif
    #endfor 
         
    [success, inputVolumeNode] = slicer.util.loadVolume( fnm, returnNode=True)

    inputFiducialNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
    inputFiducialNode.CreateDefaultDisplayNodes()
    inputFiducialNode.SetName(inputVolumeNode.GetName()+"_VertebraLocationPoints")  
    inputFiducialNode.AddFiducialFromArray(c7p)
    
    #inputFiducialNode.AddFiducial(RASpt[0],RASpt[1],RASpt[2])
    # Which vertebra
    vtID = 7
    # C7 center of mass
    # run the segmentation
    self.logic.run(inputVolumeNode, inputFiducialNode, vtID, True)    

    self.etm=time.time()
    tm=self.etm - self.stm
    print("Time: "+str(tm)+"  seconds")
    
    self.delayDisplay('Test passed!')
  #enddef

    
