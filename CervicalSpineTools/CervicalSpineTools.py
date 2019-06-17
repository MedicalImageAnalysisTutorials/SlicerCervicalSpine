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
#-------------------------------------------------------------------------------------#
#  Slicer 4.10                                                                        #    
#  Updated: 17.6.2019                                                                 #    
#=====================================================================================#

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
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing.dummy import Process  

import VisSimCommon

#TODO:
# parallel process support, check python 3 
# extend to whole spine

#===================================================================
#                           Main Class
#===================================================================

class CervicalSpineTools(ScriptedLoadableModule):
  def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        parent.title = "Cervical Spine Tools"
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
class CervicalSpineToolsWidget(ScriptedLoadableModuleWidget):
  def setup(self):
    print("=======================================================")   
    print("   VisSim Cervical Spine Tools               ")
    print("=======================================================")           
        
    ScriptedLoadableModuleWidget.setup(self)
    
    # to access logic class functions and setup global variables
    self.logic = CervicalSpineToolsLogic()
    self.vsc = VisSimCommon.VisSimCommonLogic()    
    #------------------------------------------------------------
    #                     Create the GUI interface
    #------------------------------------------------------------
    # Create main collapsible Button 
    self.mainCollapsibleBtn = ctk.ctkCollapsibleButton()
    self.mainCollapsibleBtn.setStyleSheet("ctkCollapsibleButton { background-color: DarkSeaGreen  }")
    self.mainCollapsibleBtn.text = "VisSim Cervical Spine Tools"
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
      if not hasattr(self.vsc, 'vtVars'):
         self.vsc.setGlobalVariables(1)
      #end 
      self.inputFiducialNode  = self.vsc.locateItem(self.inputVolumeNode,self.inputPointEdt, 0, self.vtID)    
      #self.onInputFiducialBtnClick()
  #enddef

  def onInputPointEdtChanged(self,point):
      self.vtID = self.vtIDCoBx.currentIndex + 1   
      self.vsc.setVtIDfromEdt(point, self.vtID)
  #enddef
    
  # extract ligaments points 
  def onLigPtsChkBxChange(self):      
      nodes = slicer.util.getNodesByClass('vtkMRMLMarkupsFiducialNode')
      ligName = "_LigPts_C" # + self.vsc.vtVars['vtID']
      self.vsc.setItemChk('ligChk',self.ligPtsChkBx.checked,ligName, nodes)

      
  #enddef
  
  def onApplyBtnClick(self):
      #check if C1,C4 and C7 points are selected
      self.runBtn.setText("...please wait")
      self.runBtn.setStyleSheet("QPushButton{ background-color: red  }")
      slicer.app.processEvents()
      self.stm=time.time()
      print("time:" + str(self.stm))
      self.timeLbl.setText("                 Time: 00:00")
    
      vtID = self.vtIDCoBx.currentIndex + 1
      inputVolumeNode = self.inputSelectorCoBx.currentNode()
      pointSelected = self.inputPointEdt.text =="[0,0,0]"
      try: 
         if (not inputVolumeNode is None) and (not pointSelected) and (not self.inputFiducialNode is None):
            # create an option to use IJK point or fidicual node
            # inputImage, FiducialPoint, vertebraID, isExteranl ,C7Editbox
            self.logic.run( inputVolumeNode ,self.inputFiducialNode )
         else:
            print("C1,C4 and C7 points are not selected !")   
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
      self.vsc.openResultsFolder()
  #enddef
#===================================================================
#                           Logic
#===================================================================
class CervicalSpineToolsLogic(ScriptedLoadableModuleLogic):
 
  def getAllVertebraePoints(self,vtIDsLst,inputFiducialNode):               
      #TODO: test on all images, the direction should be involved.
      print("Spine: Compute vertebra locations !")
      # compute if C1, C4, C7 are avaialble 
      if ((vtIDsLst[0][0]!=0) and (vtIDsLst[1][0]!=0) and(vtIDsLst[3][0]!=0) and (vtIDsLst[6][0]!=0)) : 
         # we need: 1,4 and 7
         # compute displacement between C4 and C1
         v1x= (vtIDsLst[3][0]-vtIDsLst[0][0])/4 ; print(v1x)
         v1y= (vtIDsLst[3][1]-vtIDsLst[0][1])/4 ; print(v1y)
         v1z= (vtIDsLst[3][2]-vtIDsLst[0][2])/4 ; print(v1z)
         
         # C3 location
         vtIDsLst[2]=[vtIDsLst[3][0]-v1x , vtIDsLst[3][1]-v1y , vtIDsLst[3][2]-v1z]
         inputFiducialNode.AddFiducialFromArray(vtIDsLst[2])
         inputFiducialNode.SetNthFiducialLabel(4, "C3")

         ## C2 location
         #vtIDsLst[1]=[vtIDsLst[2][0]-v1x , vtIDsLst[2][1]-v1y , vtIDsLst[2][2]-v1z]
         #self.inputFiducialNode.AddFiducialFromArray(vtIDsLst[1])
         #self.inputFiducialNode.SetNthFiducialLabel(4, "C2") # note that 4 is the index of last one added so far

         # compute displacement between C7 and C4
         v2x= (vtIDsLst[6][0]-vtIDsLst[3][0])/4
         v2y= (vtIDsLst[6][1]-vtIDsLst[3][1])/4
         v2z= (vtIDsLst[6][2]-vtIDsLst[3][2])/4
         
         # C5 location
         vtIDsLst[4]=[vtIDsLst[3][0]+v1x , vtIDsLst[3][1]+v1y , vtIDsLst[3][2]+v1z]
         inputFiducialNode.AddFiducialFromArray(vtIDsLst[4])
         inputFiducialNode.SetNthFiducialLabel(5, "C5")
         
         # C6 location
         vtIDsLst[5]=[vtIDsLst[4][0]+v1x , vtIDsLst[4][1]+v1y , vtIDsLst[4][2]+v1z]
         inputFiducialNode.AddFiducialFromArray(vtIDsLst[5])
         inputFiducialNode.SetNthFiducialLabel(6, "C6")

         fnm = os.path.join(self.vsc.vtVars['outputPath'], inputFiducialNode.GetName()+".fcsv")                             
         sR = slicer.util.saveNode(inputFiducialNode, fnm ) 
      else:
         print("Error, at least locate 4 points: C1, C2, C4, C7!!!")
         return -1
      #endif
      return vtIDsLst
  #enddef
         

  def run( self, inputVolumeNode, inputFiducialNode):       
      self.vsc   = VisSimCommon.VisSimCommonLogic()
      self.vsc.setGlobalVariables(1)
      self.inputVolumeNode   = inputVolumeNode
      self.inputFiducialNode = inputFiducialNode
      self.tableName =  inputVolumeNode.GetName()+"_tbl"
      self.outputPaths      = []
      self.intputCropPaths  = []
      self.inputPoints      = []
      self.modelPaths       = []
      vtIDsLst              = []
      missingVt             = range(1,8)
      rng= 7 # number of vertebra
      for j in range(7):
          vtIDsLst.append([0,0,0])
          self.intputCropPaths.append("")
          self.inputPoints.append([0,0,0])
      #endfor    
      for j in range(inputFiducialNode.GetNumberOfFiducials()):
           l= inputFiducialNode.GetNthFiducialLabel(j)
           k = int(l[1])
           missingVt[k-1]= 0 
           inputFiducialNode.GetNthFiducialPosition(j,vtIDsLst[k-1])
      #endfor       
               
      if sum(missingVt)>0:         
          vtIDsLst =  self.getAllVertebraePoints(vtIDsLst,inputFiducialNode)
      #endif       
           
       
      # create output paths
      for i in range(rng):
          #get other points 
          #create a list from the points
          vtID = i+1
          outputPath = self.vsc.vtVars['outputPath']+","+inputVolumeNode.GetName()+"_C"+str(vtID)
          outputPath = os.path.join(*outputPath.split(","))
          self.outputPaths.append( outputPath) 
          modelCropPath       = self.vsc.vtVars['modelPath']+ ',Default' +",Mdl" + self.vsc.vtVars['Styp']+ str(vtID)        +self.vsc.vtVars['imgType'] 
          modelCropPath       = os.path.join(*modelCropPath.split(","))
          self.modelPaths.append(modelCropPath)  
          if not os.path.exists(self.outputPaths[i]):
              os.makedirs(self.outputPaths[i])
          #endif
      #endfor  

      # try parallel but slicer crash
      # try if we collect information then run popen
      """
      i=range(rng)
      pool = ThreadPool(9)
      b=pool.map(self.runCroppingAll,i)
      pool.close()
      pool.join()
      """
      #---------------------  process ------------------------------------           
      for i in range (rng):
          vtID = i+1 
          segNodeName       = self.inputVolumeNode.GetName() + "_C"   +str(vtID)+ "S.seg"                 
          ligPtsNodeName    = self.inputVolumeNode.GetName() + "_C"   +str(vtID)+ "_LigPts"
          transNodeName     = self.inputVolumeNode.GetName() + "_C"   +str(vtID)+ "_Transform"
          self.runCroppingAll(i)
          self.runElastixAll(i)

          self.resTransPath  = os.path.join(self.outputPaths[i] ,"TransformParameters.0.txt")
          resOldDefPath = os.path.join(self.outputPaths[i] , "deformationField"+self.vsc.vtVars['imgType'])
          resDefPath    = os.path.join(self.outputPaths[i] , self.inputVolumeNode.GetName()+"_C"+str(vtID)+"_dFld"+self.vsc.vtVars['imgType'])
          #remove old result files:
          if os.path.isfile(resOldDefPath):
             os.remove(resOldDefPath) 
          if os.path.isfile(resDefPath):
             os.remove(resDefPath)
          
          self.runTransformixAll(i)

          os.rename(resOldDefPath,resDefPath)

          [success, vtTransformNode] = slicer.util.loadTransform(resDefPath, returnNode = True)
          vtTransformNode.SetName(transNodeName)

           
          modelSegPath   = self.vsc.vtVars['modelPath'] + self.vsc.vtVars['segT']+",Mdl"+self.vsc.vtVars['Styp']+ str(vtID)+ self.vsc.vtVars['sgT'] +self.vsc.vtVars['imgType'] 
          modelSegPath   =   os.path.join(*modelSegPath.split(","))
           
          [success, vtResultSegNode] = slicer.util.loadSegmentation(modelSegPath, returnNode = True)
          vtResultSegNode.SetName(segNodeName)
          vtResultSegNode.SetAndObserveTransformNodeID(vtTransformNode.GetID()) # movingAllMarkupNode should be loaded, the file contains all points
          slicer.vtkSlicerTransformLogic().hardenTransform(vtResultSegNode) 
          vtResultSegNode.CreateClosedSurfaceRepresentation() 
          fnm = os.path.join(self.outputPaths[i] , vtResultSegNode.GetName()+".nrrd")                             
          sR = slicer.util.saveNode(vtResultSegNode, fnm ) 

          if self.vsc.s2b(self.vsc.vtVars['ligChk']):            
             modelCropImgLigPtsPath = self.vsc.vtVars['modelPath'] +self.vsc.vtVars['vtPtsLigDir']+","+self.vsc.vtVars['Styp']+ str(vtID)+self.vsc.vtVars['vtPtsLigSuff']+".fcsv"
             modelCropImgLigPtsPath        = os.path.join(*modelCropImgLigPtsPath.split(","))
             [success, vtLigPtsNode] = slicer.util.loadMarkupsFiducialList  (modelCropImgLigPtsPath, returnNode = True)
             vtLigPtsNode.GetDisplayNode().SetSelectedColor(1,0,0)           
             vtLigPtsNode.GetDisplayNode().SetTextScale(0.5)
             vtLigPtsNode.SetName(ligPtsNodeName)
             vtLigPtsNode.SetAndObserveTransformNodeID(vtTransformNode.GetID()) # movingAllMarkupNode should be loaded, the file contains all points
             slicer.vtkSlicerTransformLogic().hardenTransform(vtLigPtsNode) # apply the transform
             # needed in extract scaled model
             self.vtLigPtsNode = vtLigPtsNode
             fnm = os.path.join(self.outputPaths[i] , vtLigPtsNode.GetName()+".fcsv")                             
             sR = slicer.util.saveNode(vtLigPtsNode, fnm ) 
          #endif 

          self.masterNode = slicer.util.getNode(self.inputVolumeNode.GetName() + "_C"+str(vtID)         )
          self.segNode    = slicer.util.getNode(segNodeName  )
          print(self.tableName)
          print( self.masterNode.GetName()[0:-3]+"_tbl")
          try:
             # if table exists, don't create a new one.
             self.spTblNode =  slicer.util.getNode(self.tableName)
             print("table is found .........................................")
          except:
             #table will be created in getItemInfo                
             self.spTblNode= None
             print("table is NOT found .........................................")
          #endtry
          #return -1
          self.getVertebraInfoAll( i )
      #endfor

      self.spTblNode.RemoveRow(0)   
     
      self.vsc.vtVars['segNodeCoM']= self.vsc.vtVars['segNodeCoM']
      self.vtResultSegNode = vtResultSegNode
      #Remove temporary files and nodes:
      self.vsc.removeTmpsFiles()    
      return vtResultSegNode
  #enddef


  def runCroppingAll(self,i):
      vtID = i+1
      print("------------ cropping process: " +str(i))
      # try crop in parallel:   
      for j in range (self.inputFiducialNode.GetNumberOfFiducials() ):
          if self.inputFiducialNode.GetNthFiducialLabel(j)==("C"+str(vtID)):
             break
          #endif
      #endfor
      self.inputPoints[i]=self.vsc.ptRAS2IJK(self.inputFiducialNode,self.inputVolumeNode,j)
      print(self.inputPoints[i])
      self.inputPoints[i] ="["+str(self.inputPoints[i][0])+","+str(self.inputPoints[i][1])+","+str(self.inputPoints[i][2])+"]"
      self.intputCropPaths[i]=self.vsc.runCropping(self.inputVolumeNode, self.inputPoints[i], self.vsc.vtVars['croppingLength'], self.vsc.vtVars['RSxyz'],self.vsc.vtVars['hrChk'], str(  vtID ) )
      
      print("ii = " +str(i)+ "  img :"+ self.intputCropPaths[i]          )
      print("------------ cropping process: " +str(i)+"  ...complete")
  #enddef

  def runElastixAll(self,i):
      vtID = i+1
      self.vsc.runElastix(self.vsc.vtVars['elastixBinPath'], self.intputCropPaths[i], self.modelPaths[i], self.outputPaths[i], self.vsc.vtVars['parsPath'] ,self.vsc.vtVars['noOutput'], "674")
  #enddef

  def runTransformixAll(self,i):
      vtID = i+1
      self.vsc.runTransformix(self.vsc.vtVars['transformixBinPath'] ,self.modelPaths[i], self.outputPaths[i], self.resTransPath, self.vsc.vtVars['noOutput'], "1254")
  #enddef
                
  def getVertebraInfoAll(self,i):
      vtID = i+1
      self.spTblNode = self.vsc.getItemInfo( self.segNode  , self.masterNode , self.spTblNode, vtID)      
      fnm = os.path.join(self.vsc.vtVars['outputPath'], self.spTblNode.GetName()+".tsv")                             
      sR = slicer.util.saveNode(self.spTblNode, fnm ) 
  #enddef
                              
#===================================================================
#                           Test
#===================================================================
class CervicalSpineToolsTest(ScriptedLoadableModuleTest):

  def setUp(self):
      slicer.mrmlScene.Clear(0)
      self.logic =  CervicalSpineToolsLogic()
      self.vsc   =  VisSimCommon.VisSimCommonLogic()
  #endef      

  def runTest(self):
      self.setUp()
      self.testSlicerCervicalSpineTools()
  #endef      

  def testSlicerCervicalSpineTools(self):
      self.delayDisplay("Starting testSlicerCervicalSpineTools test")
      self.stm=time.time()

      self.vsc.vtVars = self.vsc.setGlobalVariables(1)

      tstImg = 1 # 1 =CT, 2 = MR 
      if tstImg ==1:
         imgWebLink = "https://cloud.uni-koblenz-landau.de/s/Mb6JHLdWw5MEPB2/download"
         fnm = os.path.join(*(self.vsc.vtVars['vissimPath'] +",imgCT"+self.vsc.vtVars['imgType']).split(","))
         #  CT: define a markup with all locations  
         c1p = [-5.507 ,  -20.202 , -18.365 ]
         c2p = [-1.169  , -20.089 , -45.021 ]
         c4p = [ 0.390  , -11.754 , -74.194 ]
         c7p = [ 0.390  , -26.445 ,-115.821 ]
      else:
         imgWebLink = "https://cloud.uni-koblenz-landau.de/s/ieyDfHpCjHNpZXi/download"
         fnm = os.path.join(*(self.vsc.vtVars['vissimPath'] +",imgMR"+self.vsc.vtVars['imgType']).split(","))
         # MR: define a markup with all locations     
         c1p = [  1.062  ,  53.364 ,  145.951 ]
         c2p = [ -0.408  ,  44.283 ,  122.994 ]
         c4p = [ -0.408  ,  36.107 ,   89.941 ]
         c7p = [ -0.408  ,   3.439 ,   54.432 ]
      #endif
        
      # don't download if already downloaded      
      print(fnm)                 
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
        
      [success, inputVolumeNode] = slicer.util.loadVolume( fnm, returnNode=True)
    
      inputFiducialNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsFiducialNode")
      inputFiducialNode.CreateDefaultDisplayNodes()
      #TODO: change the name to points 
      inputFiducialNode.SetName("VertebraLocationPoint")  
      inputFiducialNode.AddFiducialFromArray(c1p)
      inputFiducialNode.SetNthFiducialLabel(0, "C1")
      inputFiducialNode.AddFiducialFromArray(c2p)
      inputFiducialNode.SetNthFiducialLabel(1, "C2")
      inputFiducialNode.AddFiducialFromArray(c4p)
      inputFiducialNode.SetNthFiducialLabel(2, "C4")
      inputFiducialNode.AddFiducialFromArray(c7p)
      inputFiducialNode.SetNthFiducialLabel(3, "C7")

      vtSegNode = self.logic.run(inputVolumeNode, inputFiducialNode )
    
      # Display
      self.vsc.dispSeg(inputVolumeNode,vtSegNode,34) # 34: 4up table layout

      self.etm=time.time()
      tm=self.etm - self.stm
      print("Time: "+str(tm)+"  seconds")

      self.delayDisplay('Test testSlicerCervicalSpineTools passed!')
  #enddef

 
