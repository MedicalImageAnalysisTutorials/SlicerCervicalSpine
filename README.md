**Cervical Spine Image Analysis**

<img src="https://github.com/MedicalImageAnalysisTutorials/SlicerCervicalSpine/blob/master/CervicalSpine.png" width="400" height="400">

This is a [3D Slicer](https://github.com/Slicer/Slicer) plugin that uses [elastix toolbox](https://github.com/SuperElastix/elastix) for Multi-modal CervicalSpine Images segmentation and analysis. More information can be found [here](https://mtixnat.uni-koblenz.de). The elastix parameters file can be found [here](http://elastix.bigr.nl/wiki/index.php/Par0053)

**Tested on:** 
Slicer 4.10.1, Windows 10 and Ubuntu 18.04 

This project contains two modules:

  1. Cervical Vertebra Tool.
  2. Cervical Spine Tool. 

The first tool works on a single vertebra. the second tool works on all vertebrae. Both tools do segmentation, origin and insertion points detection of the cervical spine vertebrae ligaments.

These tools can be easily extended to solve similar problems or adding more features. Note that the segmentation accuracy is not high, hence it can be used as an estimated segmentation or as an input to another enhanced algorithm. 

Please cite our papers:
*  Ibraheem AL-Dhamari, Sabine Bauer, Eva Keller and Dietrich Paulus (2019), Automatic Detection of Cervical Spine Ligaments Origin and Insertion Points. Accepted in The IEEE International Symposium on Biomedical Imaging (ISBI) 2019, Venice, Italy.

*  Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, (2018), Automatic Multi-modal Cervical Spine Image Atlas Segmentation Using Adaptive Stochastic Gradient Descent, Bildverarbeitung für die Medizin 2018 pp 303-308.

Please share your cervical spine dataset with us. 

Your contribution is welcome! 




