REM Create links to github sources:
Rem /home/%USERNAME%/.config/NA-MIC/Extensions-28257
REM open cmd as admin then run: mklink dest src
set ext="SlicerCervicalSpine"
set vr="\lib\Slicer-4.10"
set p1="C:\Users\%USERNAME%\AppData\Roaming\NA-MIC\Extensions-28257\"
set p2="\qt-scripted-modules\"
set src="D:\myGitHub\"%ext%"\4.10\"

rename %p1%%ext%%vr%%p2%"CervicalVertebraTools.py" "CervicalVertebraTools.py.bk"
mklink %p1%%ext%%vr%%p2%"CervicalVertebraTools.py" %src%"CervicalVertebraTools""\CervicalVertebraTools.py"

rename %p1%%ext%%vr%%p2%"CervicalSpineTools.py" "CervicalSpineTools.py.bk"
mklink %p1%%ext%%vr%%p2%"CervicalSpineTools.py" %src%"CervicalSpineTools""\CervicalSpineTools.py"

