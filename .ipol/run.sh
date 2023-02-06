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
PLOTFILE=$6
    
echo "INPUT3DDom = $INPUT3DDom"

if test -f "$INPUT"
then
  echo "----------------------------------------"
  echo "-----Generating 2D ---------------------"
  echo "----------------------------------------"

  convert ${INPUT} input.pgm
  ${EXEC} -n ${NBTERM} -a ${APERF}  -d input.pgm -x graphExport.xml >> algo_info.txt
  set $(identify -format '%w %h' ${INPUT})
  width=$1
  height=$2
  convert -density 400 -resize ${width}x${height}  -crop ${width}x${height} result.svg result.png
  echo "algoDim=2" >> algo_info.txt 
fi

if test -f "$INPUT3D"
then 
  echo "----------------------------------------"
  echo "-----Generating 3D ---------------------"
  echo "----------------------------------------"

  ${EXEC3D} -n ${NBTERM} -a ${APERF}  -o result.obj -x graphExport.xml
  key=$(basename $(pwd))
  demo_id=$(basename $(dirname $(pwd)))
  viewer_url="https://3dviewer.net#https://ipolcore.ipol.im/api/core/shared_folder/run/${demo_id}/${key}/result.obj,https://ipolcore.ipol.im/api/core/shared_folder/run/${demo_id}/${key}/result.mtl"
  iframe="<iframe id='3dviewerplayer' type='text/html' width='620' height='460' src='$viewer_url' "
  iframe="$iframe frameborder='5' scrolling='no' allowfullscreen webkitallowfullscreen mozallowfullscreen></iframe>"
  echo "url=$iframe" >> algo_info.txt
  echo "algoDim=3" >> algo_info.txt 
fi




if test -f "$INPUT3DDom"
then
  echo "----------------------------------------"
  echo "-----Generating 3D  Dom-----------------"
  echo "----------------------------------------"
  ${EXEC3D} -n ${NBTERM} -a ${APERF} -d $INPUT3DDom -o resultVessel.obj -x graphExport.xml
  volBoundary2obj $INPUT3DDom liver05Domain.obj
  mergeObj resultVessel.obj liver05Domain.obj result.obj --nameGrp1  vessel --nameGrp2  liver  --materialOne 0.7 0.2 0.2 1.0 --materialTwo +0.4  0.4 0.5 0.2 
  key=$(basename $(pwd))
  demo_id=$(basename $(dirname $(pwd)))
  viewer_url="https://3dviewer.net#https://ipolcore.ipol.im/api/core/shared_folder/run/${demo_id}/${key}/result.obj,https://ipolcore.ipol.im/api/core/shared_folder/run/${demo_id}/${key}/result.mtl"
  iframe="<iframe id='3dviewerplayer' type='text/html' width='620' height='460' src='$viewer_url' "
  iframe="$iframe frameborder='5' scrolling='no' allowfullscreen webkitallowfullscreen mozallowfullscreen></iframe>"
  echo "url=$iframe" >> algo_info.txt
  echo "algoDim=3" >> algo_info.txt 
fi


  echo "----------------------------------------"
  echo "-----Generating Stat Radius Curves -----"
  echo "----------------------------------------"
  COMMANDStat1="xml2graph graphExport.xml"
  COMMANDStat2="graph2statBifRad vertex.txt edges.txt radius.txt stat.dat"
  COMMANDStat3="gnuplot $PLOTFILE"

  for c in COMMANDStat1 COMMANDStat2 COMMANDStat3
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

  
  


