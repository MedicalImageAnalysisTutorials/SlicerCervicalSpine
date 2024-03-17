**Cervical Spine Image Analysis**

<img src="https://github.com/MedicalImageAnalysisTutorials/SlicerCervicalSpine/blob/master/CervicalSpine.png" width="400" height="400">

This is a [3D Slicer](https://github.com/Slicer/Slicer) plugin that uses [elastix toolbox](https://github.com/SuperElastix/elastix) for Multi-modal CervicalSpine Images segmentation and analysis. More information can be found in the related publications. The elastix parameters file can be found [here](https://github.com/MedicalImageAnalysisTutorials/SlicerCervicalSpine/tree/master/pars)

**Tested on:** 
Slicer 5.6, Windows 10 and Ubuntu 20.04 

This project contains two modules:

  1. Cervical Vertebra Tool.
  2. Cervical Spine Tool. 

The first tool works on a single vertebra. the second tool works on all vertebrae. Both tools do segmentation, origin and insertion points detection of the cervical spine vertebrae ligaments.

These tools can be easily extended to solve similar problems or adding more features. Note that the segmentation accuracy is not high, hence it can be used as an estimated segmentation or as an input to another enhanced algorithm. 

Please cite our papers:
*  Ibraheem AL-Dhamari, Sabine Bauer, Eva Keller and Dietrich Paulus (2019), Automatic Detection of Cervical Spine Ligaments Origin and Insertion Points. Accepted in The IEEE International Symposium on Biomedical Imaging (ISBI) 2019, Venice, Italy.

*  Ibraheem Al-Dhamari, Sabine Bauer, Dietrich Paulus, (2018), Automatic Multi-modal Cervical Spine Image Atlas Segmentation Using Adaptive Stochastic Gradient Descent, Bildverarbeitung für die Medizin 2018 pp 303-308.


Your contribution is welcome! 

# Notes:  

* For general 3D Slicer questions, please use Slicer [forum](https://discourse.slicer.org), many experts will be able to help you there. 
* For Slicer Spine related questions, comemnts, or feedback, please use github [discussion](https://github.com/MedicalImageAnalysisTutorials/SlicerCervicalSpine/discussions/categories/q-a) section. 
* For  Slicer Spine related bugs, please open a new [issue](https://github.com/MedicalImageAnalysisTutorials/SlicerCervicalSpine/issues) if needed. Please mention your operating system, slicer version, and the error message you see in python interactor.
* For sharing private information, dataset, a future project proposal or cooperation, please use the [email](ia@idhamari.com), use SlicerSpine in the subject. 

# License

Copyright 2019, [Dr. Ibraheem AL-Dhamari](https://idhamari.com).

Licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0). 




