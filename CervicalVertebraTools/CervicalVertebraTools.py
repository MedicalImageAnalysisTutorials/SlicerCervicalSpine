#=====================================================================================#
#  3D Slicer [1] plugin that uses elastix toolbox [2] Plugin for Automatic Cervical   # 
#  Spine Image Segmentation and other useful features extraction[3]                   #
#  Locate a single vertebra and get its segmentation and other features               # 
#  More info can be found at [3 and 4]                                                #
#  Sample vertebra datasets can be downloaded using Slicer Datastore module           #
#  Note: it is recommended to use deep learning based segmentation instead of         #
#         this extension                                                              #
#                                                                                     #
#  Contributers: Ibraheem Al-Dhamari,  ia@idhamari                                    #
#  [1] https://www.slicer.org                                                         #
#  [2] http://elastix.isi.uu.nl                                                       #
#  [3] Ibraheem AL-Dhamari, Sabine Bauer, Eva Keller and Dietrich Paulus, (2019),     #
#      Automatic Detection Of Cervical Spine Ligaments Origin And Insertion Points,   #
#      IEEE International Symposium on Biomedical Imaging (ISBI), Venice, Italy.      #
#  [4] Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, (2018), Automatic          #
#      Multi-modal Cervical Spine Image Atlas Segmentation Using Adaptive Stochastic  #
#      Gradient Descent, Bildverarbeitung feur die Medizin 2018 pp 303-308.           #  
#                                                                                     #
#-------------------------------------------------------------------------------------#
#  Slicer 5.6                                                                         #    
#  Updated: 17.3.2024                                                                 #    
#=====================================================================================#

import os, re , datetime, time ,shutil, unittest, logging, zipfile, urllib.request, stat,  inspect, glob
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
import SampleData

#Dependant Modules
import VisSimCommon #Dependant Modules

# Terminology
#  img         : ITK image 
#  imgNode     : Slicer Node
#  imgName     :  Filename without the path and without extension
#  imgPath     : wholePath + Filename and extension

#===================================================================
#                           Main Class
#===================================================================
class CervicalVertebraTools(ScriptedLoadableModule):
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
    
#===================================================================
#                           Main Widget
#===================================================================
class CervicalVertebraToolsWidget(ScriptedLoadableModuleWidget):
  def setup(self):
    print("=======================================================")   
    print("   VisSim Cervical Vertebra Tools               ")
    print("=======================================================")           
    ScriptedLoadableModuleWidget.setup(self)    
    
    # to access logic class functions and setup global variables
    self.logic = CervicalVertebraToolsLogic()
    self.vsc = VisSimCommon.VisSimCommonLogic()   
    #----------------------------------------------------------------- 
    #                     Create the GUI interface
    #-----------------------------------------------------------------    
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
    # Create a textbox for vertebra location
    self.inputPointEdt = qt.QLineEdit()
    self.inputPointEdt.setFixedHeight(20)
    self.inputPointEdt.setFixedWidth(100)
    self.inputPointEdt.setText("[0,0,0]")
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

  #------------------------------------------------------------------------
  #                        Define GUI Elements Functions
  #------------------------------------------------------------------------

  # Locate a vertebra
  def onVtIDCoBxChange(self):
      self.inputVolumeNode=self.inputSelectorCoBx.currentNode()
      self.vtID = self.vtIDCoBx.currentIndex + 1   
      self.inputVolumeNode  =  self.inputVolumeNode 
      self.inputFiducialNode = None
      nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
      for f in nodes:
          if ((f.GetName() == self.inputVolumeNode.GetName()+"_vtLocations") ):
             #replace  current 
             #print("inputFiducialNode exist")
             self.inputFiducialNode = f  
             newNode= False

      if not hasattr(self.vsc, 'vtVars'):
         self.vsc.setGlobalVariables(1)

      self.inputFiducialNode =  self.vsc.locateItem(self.inputVolumeNode,self.inputPointEdt, 0, self.vtID)  
      
      #self.onInputFiducialBtnClick()

  def onInputPointEdtChanged(self,point):
      if not hasattr(self.vsc, 'vtVars'):
         self.vsc.setGlobalVariables(1)

      self.vtID = self.vtIDCoBx.currentIndex + 1   
      self.vsc.setVtIDfromEdt(point,self.vtID)

  # extract ligaments points 
  def onLigPtsChkBxChange(self):      
      nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
      ligName = "_LigPts" # + self.vsc.vtVars['vtID']
      self.vsc.setItemChk('ligChk',self.ligPtsChkBx.checked,ligName, nodes)
  
  def onApplyBtnClick(self):      
      self.runBtn.setText("...please wait")
      self.runBtn.setStyleSheet("QPushButton{ background-color: red  }")
      slicer.app.processEvents()
      self.stm=time.time()
      print("time:" + str(self.stm))
      self.timeLbl.setText("                 Time: 00:00")
      methodID =0
      vtID = self.vtIDCoBx.currentIndex + 1
      inputVolumeNode = self.inputSelectorCoBx.currentNode()
      pointSelected = self.inputPointEdt.text =="[0,0,0]"
      try:         
         if (not inputVolumeNode is None) and (not pointSelected) and (not self.inputFiducialNode is None):
            # create an option to use IJK point or fidicual node
            # inputImage, FiducialPoint, vertebraID, isExteranl ,C7Editbox
            vtSegNode =self.logic.run(inputVolumeNode ,self.inputFiducialNode, vtID, methodID )
            self.vsc.dispSeg(inputVolumeNode ,vtSegNode, 34 )
         else:
            print("error in input, please check input image or input point")
         
      except Exception as e:
            print("STOPPED: error in input or during the segmentation")
            print(e)
            
            
      self.etm=time.time()
      tm=self.etm - self.stm
      self.timeLbl.setText("Time: "+str(tm)+"  seconds")
      self.runBtn.setText("Run")
      self.runBtn.setStyleSheet("QPushButton{ background-color: DarkSeaGreen  }")
      slicer.app.processEvents()
       
  def onOpenResultFolderBtnClick(self):
      self.vsc.openResultsFolder()
    

#===================================================================
#                           Logic
#===================================================================
class CervicalVertebraToolsLogic(ScriptedLoadableModuleLogic):

   
  #--------------------------------------------------------------------------------------------        
  #                       Segmentation Process 
  #--------------------------------------------------------------------------------------------        
  # This method perform the atlas segementation steps
  def run(self, inputVolumeNode, inputFiducialNode, vtIDT, vtMethodID):
      logging.info('Processing started')
      self.vsc   = VisSimCommon.VisSimCommonLogic()
      self.vsc.setGlobalVariables(1)
      print("setup local paths variables.........................")  
      vtID =int(vtIDT) 
      # set the correct models paths:
      modelPath       =   self.vsc.vtVars['modelPath']+  ',Default'+",Mdl" + self.vsc.vtVars['Styp']+ str(vtID)        +self.vsc.vtVars['imgType'] 
      modelPath       =   os.path.join(*modelPath.split(","))
      modelSegPath    =   self.vsc.vtVars['modelPath'] + self.vsc.vtVars['segT']+",Mdl"+self.vsc.vtVars['Styp']+ str(vtID)+ self.vsc.vtVars['sgT'] +self.vsc.vtVars['imgType'] 
      modelSegPath    =   os.path.join(*modelSegPath.split(","))
      modelLigPtsPath = self.vsc.vtVars['modelLigPtsPath'] +str(vtID)+self.vsc.vtVars['vtPtsLigSuff']+".fcsv"

      # set the results paths:       
      resTransPathOld  = os.path.join(self.vsc.vtVars['outputPath'] ,"TransformParameters.0.txt")
      resTransPath=resTransPathOld[0:-6]+'Pars.txt'        
      resOldDefPath = os.path.join(self.vsc.vtVars['outputPath'] , "deformationField"+self.vsc.vtVars['imgType'])
      resDefPath    = os.path.join(self.vsc.vtVars['outputPath'] , inputVolumeNode.GetName()+"_C"+str(vtID)+"_dFld"+self.vsc.vtVars['imgType'])
      inputImgName  = inputVolumeNode.GetStorageNode().GetFileName()
      inputImgName  = basename(os.path.splitext(inputImgName)[0])    
        
      segNodeName       = inputVolumeNode.GetName() + "_C"   +str(vtID)+".Seg"                 
      ligPtsNodeName    = inputVolumeNode.GetName() + "_C"   +str(vtID)+"_LigPts"
      transNodeName     = inputVolumeNode.GetName() + "_C"   +str(vtID)+ "_Transform"

      print("removeOtputsFolderContents ........................ ")  
      #only for this vertebra
      self.vsc.removeOtputsFolderContents()
      # check if the model is found
      if not isfile(modelPath): 
            print >> sys.stderr, "ERROR: model is not found"            
            print("modelPath: " + modelPath)
            return -1
      # endif

      print("ptRAS2IJK ........................ ")          
      # for multiple vertebra, it is enough to have one markup node
      # Get IJK point from the fiducial to use in cropping
      newFid= True
      for j in range (inputFiducialNode.GetNumberOfControlPoints() ):
             if inputFiducialNode.GetNthControlPointLabel(j)==("C"+str(vtID)) :
                break

      # Get IJK point from the fiducial to use in cropping          
      inputPoint = self.vsc.ptRAS2IJK(inputFiducialNode,inputVolumeNode,j)
      # TODO: add better condition
      if  np.sum(inputPoint)== 0 :
            print("Error: select vertebra point")
            return -1

      fnm = os.path.join(self.vsc.vtVars['outputPath'] , inputVolumeNode.GetName()+"_C"+str(vtID)+"_vtLocations.fcsv")                           
      sR = slicer.util.saveNode(inputFiducialNode, fnm )  

      #Remove old resulted nodes
      for node in slicer.util.getNodes():
          if ( ligPtsNodeName    == node): slicer.mrmlScene.RemoveNode(node) 
          if ( segNodeName       == node): slicer.mrmlScene.RemoveNode(node) 
          if ( transNodeName == node): slicer.mrmlScene.RemoveNode(node) 

      inputPointT = self.vsc.v2t(inputPoint) 
      print ("************  Cropping  **********************")
      self.vsc.vtVars['intputCropPath'] = self.vsc.runCropping(inputVolumeNode, inputPointT,self.vsc.vtVars['croppingLength'],  self.vsc.vtVars['RSxyz'],  self.vsc.vtVars['hrChk'], str(vtID) )                    
      croppedNode = slicer.util.loadVolume(self.vsc.vtVars['intputCropPath'])
      croppedNode.SetName(inputVolumeNode.GetName()+"_C"+str(vtID)+"Crop")        

      print ("************  Register model to cropped input image **********************")
      cTI = self.vsc.runElastix(self.vsc.vtVars['elastixBinPath'],self.vsc.vtVars['intputCropPath'],  modelPath, self.vsc.vtVars['outputPath'], self.vsc.vtVars['parsPath'], self.vsc.vtVars['noOutput'], "554")
      copyfile(resTransPathOld, resTransPath)
      #genrates deformation field 
      cTR = self.vsc.runTransformix(self.vsc.vtVars['transformixBinPath'],modelPath, self.vsc.vtVars['outputPath'], resTransPath, self.vsc.vtVars['noOutput'], "556")
      # rename fthe file:
      os.rename(resOldDefPath,resDefPath)       
      print ("************  Load deformation field Transform  **********************")
      vtTransformNode = slicer.util.loadTransform(resDefPath)
      vtTransformNode.SetName(transNodeName)
      print ("************  Transform segmentation  **********************")
      vtSegNode = slicer.util.loadSegmentation(modelSegPath)
      vtSegNode.SetName(segNodeName)
      vtSegNode.SetAndObserveTransformNodeID(vtTransformNode.GetID()) # movingAllMarkupNode should be loaded, the file contains all points
      slicer.vtkSlicerTransformLogic().hardenTransform(vtSegNode) # apply the transform
      vtSegNode.CreateClosedSurfaceRepresentation() 
      fnm = os.path.join(self.vsc.vtVars['outputPath'] , vtSegNode.GetName()+".nrrd")                             
      sR = slicer.util.saveNode(vtSegNode, fnm )  
                     
      if self.vsc.s2b(self.vsc.vtVars['ligChk']):            
           print ("************  Transform Ligaments Points **********************")
           #vtLigPtsNode = slicer.util.loadMarkupsFiducialList  (modelLigPtsPath)
           vtLigPtsNode = slicer.util.loadMarkups (modelLigPtsPath)
           vtLigPtsNode.GetDisplayNode().SetSelectedColor(1,0,0)           
           vtLigPtsNode.GetDisplayNode().SetTextScale(0.5)                      
           vtLigPtsNode.SetName(ligPtsNodeName)
           vtLigPtsNode.SetAndObserveTransformNodeID(vtTransformNode.GetID()) 
           slicer.vtkSlicerTransformLogic().hardenTransform(vtLigPtsNode) # apply the transform
           # needed in extract scaled model
           self.vtLigPtsNode = vtLigPtsNode    
           fnm = os.path.join(self.vsc.vtVars['outputPath'] , vtLigPtsNode.GetName()+".fcsv")                             
           sR = slicer.util.saveNode(vtLigPtsNode, fnm ) 
   
      # Display the result if no error
      # Clear vertebra location labels
      if  (cTI==0) and (cTR==0):
          # change the model type from vtk to stl 
          msn=slicer.vtkMRMLModelStorageNode()
          msn.SetDefaultWriteFileExtension('stl')
          slicer.mrmlScene.AddDefaultNode(msn)        
          print("get vertebra information")
          tableName =  inputVolumeNode.GetName()+"_tbl"
          # create only if it does not exist
          try:
              # if table exists, don't create a new one.
              spTblNode =  slicer.util.getNode(tableName)
              print("table is found .........................................")
          except:
              #create a new table                
              spTblNode= None

          spTblNode = self.vsc.getItemInfo( vtSegNode, croppedNode, spTblNode, vtID)
      else:
          print("error happened during segmentation ")

        
      #Remove temporary files and nodes:
      self.vsc.removeTmpsFiles()     
      print("================= vertebra analysis is complete  =====================")
      logging.info('Processing completed')
      return vtSegNode          
                                     
#===================================================================
#                           Test
#===================================================================
class CervicalVertebraToolsTest(ScriptedLoadableModuleTest):
  
  def setUp(self):
      # this function is not called from external modules  
      slicer.mrmlScene.Clear(0)    
      # communicate with logic classes of  this extension and VisSimCommon extensions. 
 
  def runTest(self):
      self.setUp()
      self.testSlicerCervicalVertebraTools()

  def testSlicerCervicalVertebraTools(self, imgPath=None , inputPoint=None , vtID=None, methodID=None):

      self.delayDisplay("Starting testSlicerCervicalVertebraTools")

      # record duration of the test    
      self.stm=time.time()

      self.vsc   = VisSimCommon.VisSimCommonLogic()   
      self.vsc.vtVars = self.vsc.setGlobalVariables(1)
      self.logic = CervicalVertebraToolsLogic()
 
      if methodID is None:
          methodID =0 # default

      # Which vertebra
      if vtID is None:
          vtID = 7

      if imgPath is None:
         tstImg = 2 # 1 =CT, 2 = MR 
         if tstImg ==1:
            #TODO: add alternative image links 
            nodeNames='Bc11702'   
            fileNames='Bc11702.nrrd'
            urisUniKo     ="https://cloud.uni-koblenz-landau.de/s/Mb6JHLdWw5MEPB2/download"
            urisGitHub   = 'https://github.com/MedicalImageAnalysisTutorials/VisSimData/raw/master/Bc11702.nrrd'
            uris = urisGitHub          
            checksums='SHA256:f2e6623cf11566179291e648982b46a3bc6aba9abe388e24fda57e54de98eb7c'
            c7pIJK   = [ 146  , 164 , 19 ]
         else:
            nodeNames='D0040100402_3D'
            fileNames='D0040100402_3D.nrrd'
            urisUniKo     ="https://cloud.uni-koblenz-landau.de/s/ieyDfHpCjHNpZXi/download"
            urisGitHub   = 'https://github.com/MedicalImageAnalysisTutorials/VisSimData/raw/master/D0040100402_3D.nrrd'
            uris = urisGitHub          
            checksums='SHA256:a034ae045e16bdb1356e6cd9eec32ad2e3744f7f849a24abfe56cce5781dec98'
            c7pIJK = [ 253  ,   267 ,   31 ] # [-0.1703168622531166, 1.2642467778633772, 56.48611368402197]

         tmpVolumeNode =  SampleData.downloadFromURL(uris, fileNames, nodeNames, checksums )[0]
         imgPath  =  os.path.join(slicer.mrmlScene.GetCacheManager().GetRemoteCacheDirectory(),fileNames)
         slicer.mrmlScene.RemoveNode(tmpVolumeNode)
      else:
           nodeNames = os.path.splitext(os.path.basename(imgPath))[0]
 
      inputVolumeNode  = slicer.util.loadVolume(imgPath)
      inputVolumeNode.SetName(nodeNames)

      if inputPoint is None:
          inputPoint =c7pIJK
  
      #this will cause image load 
      c7p = self.vsc.ptIJK2RAS(inputPoint,imgPath)   

      # remove contents of output folder
      self.vsc.removeOtputsFolderContents()

      # create a point node for the vertebra location, this is used for cropping later
      inputFiducialNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
      inputFiducialNode.CreateDefaultDisplayNodes()
      inputFiducialNode.SetName(inputVolumeNode.GetName()+"_vtLocations")  
      inputFiducialNode. AddControlPoint(c7p)
    
      # run the segmentation
      vtSegNode = self.logic.run(inputVolumeNode, inputFiducialNode, vtID, methodID)    

      # get segmentation information
      spTblNode = self.vsc.getItemInfo(vtSegNode, inputVolumeNode, None, vtID )

      #display:
      try:
          self.vsc.dispSeg(inputVolumeNode,vtSegNode,34) # 34: 4up table layout
      except Exception as e:
             print("Can not display results! probably an external call ...")
             print(e)   

      self.etm=time.time()
      tm=self.etm - self.stm
      print("Time: "+str(tm)+"  seconds")
    
      self.delayDisplay('Test testSlicerCervicalVertebraTools passed!')