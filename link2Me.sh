# this link the installed extension source file to this repository.

ewd1=$HOME"/.config/NA-MIC/Extensions-28283"
ewd2="lib/Slicer-4.11/qt-scripted-modules"
ghb=""

#CervicalSpine
ext="SlicerCervicalSpine"
mv $ewd1/$ext/$ewd2/CervicalVertebraTools.py $ewd1/$ext/$ewd2/CervicalVertebraTools.py.bk
ln -s $PWD/CervicalVertebraTools/CervicalVertebraTools.py $ewd1/$ext/$ewd2/CervicalVertebraTools.py  
mv $ewd1/$ext/$ewd2/CervicalSpineTools.py $ewd1/$ext/$ewd2/CervicalSpineTools.py.bk
ln -s $PWD/CervicalSpineTools/CervicalSpineTools.py $ewd1/$ext/$ewd2/CervicalSpineTools.py  

