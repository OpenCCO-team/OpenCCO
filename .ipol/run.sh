#!/bin/bash
PATH=$PATH:/usr/bin/:/usr/local/bin:/opt/local/bin


NBTERM=$1
APERF=$2
INPUT=$3
INPUT3D=$4
INPUT3DDom=$5
INPUTBASE=$(basename $INPUT)
EXEC=generateTree2D
EXEC3D=generateTree3D
IPOLDIR=$6


function applyCommand
{
  for c in $*
  do
      echo "Starting command: ${!c}  "
      eval ${!c}
      if [ $? -ne 0 ] 
      then
          exit 1
      else
          echo "[done]"
      fi
  done
}


echo "INPUT3DDom = $INPUT3DDom"

if test -f "$INPUT"
then
  echo "----------------------------------------"
  echo "-----Generating 2D ---------------------"
  echo "----------------------------------------"

  COMMANDGem2D1="convert ${INPUT} input.pgm"
  COMMANDGem2D2="${EXEC} -n ${NBTERM} -a ${APERF}  -d input.pgm -x graphExport.xml >> algo_info.txt"
  set $(identify -format '%w %h' ${INPUT})
  width=$1; height=$2
  COMMANDGem2D3="convert -density 400 -resize ${width}x${height}  -crop ${width}x${height} result.svg result.png"
  applyCommand COMMANDGem2D1 COMMANDGem2D2 COMMANDGem2D3
  echo "algoDim=2" >> algo_info.txt 
elif test ! -s "$INPUT3DDom"
then 
  echo "----------------------------------------"
  echo "-----Generating 3D ---------------------"
  echo "----------------------------------------"
  COMMANDGem3D1="${EXEC3D} -n ${NBTERM} -a ${APERF}  -o result.obj -x graphExport.xml"
  applyCommand COMMANDGem3D1
  key=$(basename $(pwd))
  demo_id=$(basename $(dirname $(pwd)))
  viewer_url="https://3dviewer.net#https://ipolcore.ipol.im/api/core/shared_folder/run/${demo_id}/${key}/result.obj,https://ipolcore.ipol.im/api/core/shared_folder/run/${demo_id}/${key}/result.mtl"
  iframe="<iframe id='3dviewerplayer' type='text/html' width='620' height='460' src='$viewer_url' "
  iframe="$iframe frameborder='5' scrolling='no' allowfullscreen webkitallowfullscreen mozallowfullscreen></iframe>"
  echo "url=$iframe" >> algo_info.txt
  echo "algoDim=3" >> algo_info.txt
  echo "domain=0" >> algo_info.txt

  export DISPLAY=:1; Xvfb "$DISPLAY" -screen 0 1024x768x24 &
  meshViewer result.obj -b 240 240 240 2>&1 |  while read -r line; do if [[ $line == "[done]." ]] ; then sleep 0.5; import -window root -display :1 -screen extractVisu.png ; convert extractVisu.png  -crop 800x600+0+0  visuMeshArchive.png; pkill meshViewer;  fi; done;

  
else
  echo "----------------------------------------"
  echo "-----Generating 3D  Dom-----------------"
  echo "----------------------------------------"
  COMMANDGem3Ddom_1="${EXEC3D} -n ${NBTERM} -a ${APERF} -d $INPUT3DDom -o resultVessel.obj -x graphExport.xml"
  COMMANDGem3Ddom_2="volBoundary2obj $INPUT3DDom liver05Domain.obj"
  COMMANDGem3Ddom_3="mergeObj resultVessel.obj liver05Domain.obj result.obj --nameGrp1  vessel --nameGrp2  liver  --materialOne 0.7 0.2 0.2 1.0 --materialTwo +0.4  0.4 0.5 0.2"
  applyCommand COMMANDGem3Ddom_1 COMMANDGem3Ddom_2 COMMANDGem3Ddom_3
  key=$(basename $(pwd))
  demo_id=$(basename $(dirname $(pwd)))
  viewer_url="https://3dviewer.net#https://ipolcore.ipol.im/api/core/shared_folder/run/${demo_id}/${key}/result.obj,https://ipolcore.ipol.im/api/core/shared_folder/run/${demo_id}/${key}/result.mtl"
  iframe="<iframe id='3dviewerplayer' type='text/html' width='620' height='460' src='$viewer_url' "
  iframe="$iframe frameborder='5' scrolling='no' allowfullscreen webkitallowfullscreen mozallowfullscreen></iframe>"
  echo "url=$iframe" >> algo_info.txt
  echo "algoDim=3" >> algo_info.txt
  echo "domain=1" >> algo_info.txt
  export DISPLAY=:1; Xvfb "$DISPLAY" -screen 0 1024x768x24 &
  meshViewer resultVessel.obj -b 240 240 240  2>&1 |  while read -r line; do if [[ $line == "[done]." ]] ; then sleep 0.5; import -window root -display :1 -screen extractVisu.png ; convert extractVisu.png  -crop 800x600+0+0  visuMeshArchive.png; pkill meshViewer;  fi; done;
  meshViewer liver05Domain.obj -b 240 240 240 -c  2>&1 |  while read -r line; do if [[ $line == "[done]." ]] ; then sleep 0.5; import -window root -display :1 -screen extractVisu.png ; convert extractVisu.png  -crop 800x600+0+0  visuMeshDomain.png; pkill meshViewer;  fi; done;

fi





  echo "----------------------------------------"
  echo "-----Generating Stat Radius Curves -----"
  echo "----------------------------------------"
  COMMANDStat1="xml2graph graphExport.xml"
  COMMANDStat2="graph2statBifRad vertex.dat edges.dat radius.dat stat.dat"
  COMMANDStat3="gnuplot ${IPOLDIR}/helpers/plotStatRadius.plt"
  cp ${IPOLDIR}/helpers/readme.txt ./
  COMMANDStat4="tar cvzf graphExport.tar.gz vertex.dat edges.dat radius.dat readme.txt"
  applyCommand COMMANDStat1 COMMANDStat2 COMMANDStat3 COMMANDStat4
  

  
