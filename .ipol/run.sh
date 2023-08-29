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
FIRSTSEG=$7
X0=$8
Y0=$9
Z0=${10}
MINDISTBORDER=${11}
IMPLICITETYPE=0
IMPLICITEDIM=0
INPUTNAME=${orig_input_0}
echo "------"
echo "name : ${INPUTNAME}"
echo "------"

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

if test ! -s "$INPUT3DDom"
then
    INPUTDIM=3
    IMPLICITETYPE=0
fi

if test -f "$INPUT"
then
    INPUTDIM=2
   IMPLICITETYPE=0
fi

case $INPUTNAME in
    ##  2D shape 1 
    "e97fcd97f2d3a0c0e450e2f3c5b5eab401410d3c.png")
      IMPLICITETYPE=0
      INPUTDIM=2      
    ;;

    ##  2D shape 2
    "cd256103ba8192575e1b207045aaff03ad426e3c.png")
      IMPLICITETYPE=0
      INPUTDIM=2      
    ;;

    ## 3D liver domain    
    "3cdd37cb14ff74bb854e78ca6b9332da82a54089.png")
      IMPLICITETYPE=0
      INPUTDIM=3      
      ;;

    ## 3D toy domain    
    "ada61fd93c4a156cac2e121072c89632571d5296.png")
      IMPLICITETYPE=0
      INPUTDIM=3      
      ;;
    
    ##  3D 2 balls domain
    "31d51c4f87273fe8c0f81d938cfbe608c8a87e46.png")
      IMPLICITETYPE=0
      INPUTDIM=3      
      ;;

    ## square (implicit)
    "05704513be217481f51ba5f7ac5d8b6022c90b52.png")
      IMPLICITETYPE=2
      INPUTDIM=2      
      ;;
    ## disk (implicit)
    "babc4e2475a64382cb224402660e5a2ae2221739.png")
      IMPLICITETYPE=1
      INPUTDIM=2      
      ;;
    ## ball (implicit)
    "7c636f1832132a7bc6b737934cee8b0d0a547cc5.png")
      IMPLICITETYPE=1
      INPUTDIM=3    
      ;;
    ## box (implicit)
    "c891c555440275930d1c831acff5260a76e10d06.png")
     IMPLICITETYPE=2
     INPUTDIM=3  
     ;;
    

esac


echo "Input dim = ${INPUTDIM}"
echo "Implicit type  = ${IMPLICITETYPE}"

echo "INPUT3DDom = $INPUT3DDom"
if [ ${INPUTDIM} -eq 2 ] && [ $IMPLICITETYPE -ne 0 ]
then
  echo "----------------------------------------"
  echo "-----Generating IMPLICIT 2D ------------"
  echo "----------------------------------------"
  
  COMMANDGem2D2="${EXEC} -n ${NBTERM} -a ${APERF}  -x graphExport.xml -m ${MINDISTBORDER} >> algo_info.txt"
  if [ $IMPLICITETYPE -eq 2 ]
  then
      COMMANDGem2D2="$COMMANDGem2D2 -s"
  fi
  echo "algoImp=1" >> algo_info.txt 
  width=480; height=480;
  COMMANDGem2D3="convert -density 800 -resize ${width}x${height}  -crop ${width}x${height} result.svg result.png"
  applyCommand COMMANDGem2D1 COMMANDGem2D2 COMMANDGem2D3
  echo "algoDim=2" >> algo_info.txt 
elif test -f "$INPUT" 
then
  echo "----------------------------------------"
  echo "-----Generating 2D ---------------------"
  echo "----------------------------------------"
  echo "algoImp=0" >> algo_info.txt 

  COMMANDGem2D1="convert ${INPUT} input.pgm"
  if [ $FIRSTSEG -eq 1 ]
     then 
  COMMANDGem2D2="${EXEC} -n ${NBTERM} -a ${APERF}  -d input.pgm -x graphExport.xml -m ${MINDISTBORDER} >> algo_info.txt"
  else
  COMMANDGem2D2="${EXEC} -n ${NBTERM} -p $X0 $Y0 -a ${APERF}  -d input.pgm -x graphExport.xml -m ${MINDISTBORDER} >> algo_info.txt" 
  fi
  set $(identify -format '%w %h' ${INPUT})
  width=$1; height=$2
  COMMANDGem2D3="convert -density 400 -resize ${width}x${height}  -crop ${width}x${height} result.svg result.png"
  applyCommand COMMANDGem2D1 COMMANDGem2D2 COMMANDGem2D3
  echo "algoDim=2" >> algo_info.txt 
elif [ $INPUTDIM -eq 3 ] &&  [ $IMPLICITETYPE -ne 0 ]
then 
  echo "----------------------------------------"
  echo "-----Generating 3D IMPLICIT ---------------------"
  echo "----------------------------------------"
  COMMANDGem3D1="${EXEC3D} -n ${NBTERM} -a ${APERF}  -o result.obj -x graphExport.xml"
  if [ $IMPLICITETYPE -eq 2 ]
  then
      COMMANDGem3D1="$COMMANDGem3D1 -s"
  fi
  
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
  meshViewer result.obj -b 240 240 240 2>&1 |  while read -r line; do if [[ $line == "[done]." ]] ; then sleep 1; import -window root -display :1 -screen extractVisu.png ; convert extractVisu.png  -crop 800x600+0+0  visuMeshArchive.png; pkill meshViewer;  fi; done;

  
else
  echo "----------------------------------------"
  echo "-----Generating 3D  Dom-----------------"
  echo "----------------------------------------"
if [ $FIRSTSEG -eq 1 ]
then   
    COMMANDGem3Ddom_1="${EXEC3D} -n ${NBTERM} -a ${APERF} -m ${MINDISTBORDER} -d $INPUT3DDom -o resultVessel.obj -x graphExport.xml"
else
    COMMANDGem3Ddom_1="${EXEC3D} -n ${NBTERM} -p $X0 $Y0 $Z0 -a ${APERF} -m ${MINDISTBORDER} -d $INPUT3DDom -o resultVessel.obj -x graphExport.xml"
fi

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
  meshViewer resultVessel.obj -b 240 240 240 -l  2>&1 |  while read -r line; do if [[ $line == "[done]." ]] ; then sleep 1; import -window root -display :1 -screen extractVisu.png ; convert extractVisu.png  -crop 800x600+0+0  visuMeshArchive.png; pkill meshViewer;  fi; done;
  meshViewer liver05Domain.obj --customColorMesh 200 50 50 50 0 0 0 0 -b 240 240 240 -c -l  2>&1 |  while read -r line; do if [[ $line == "[done]." ]] ; then sleep 1; import -window root -display :1 -screen extractVisu.png ; convert extractVisu.png  -crop 800x600+0+0  visuMeshDomain.png; pkill meshViewer;  fi; done;
  meshViewer liver05Domain.obj resultVessel.obj --customColorMesh 200 50 50 50 0 0 0 0 50 50 200 255 0 0 0 0 -b 240 240 240 -c -l  2>&1 |  while read -r line; do if [[ $line == "[done]." ]] ; then sleep 2; import -window root -display :1 -screen extractVisu.png ; convert extractVisu.png  -crop 800x600+0+0  visuMeshDomainRes.png; pkill meshViewer;  fi; done;

fi





  echo "----------------------------------------"
  echo "-----Generating Stat Radius Curves -----"
  echo "----------------------------------------"
  COMMANDStat1="xml2graph graphExport.xml"
  COMMANDStat2="graph2statBifRad stat.dat vertex.dat edges.dat radius.dat"
  COMMANDStat3="gnuplot ${IPOLDIR}/helpers/plotStatRadius.plt"
  cp ${IPOLDIR}/helpers/readme.txt ./
  COMMANDStat4="tar cvzf graphExport.tar.gz vertex.dat edges.dat radius.dat readme.txt"
  applyCommand COMMANDStat1 COMMANDStat2 COMMANDStat3 COMMANDStat4
  

  
