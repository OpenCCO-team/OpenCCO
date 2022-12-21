#!/bin/bash
PATH=$PATH:/usr/bin/:/usr/local/bin:/opt/local/bin


NBTERM=$1
APERF=$2
INPUT=$3
INPUT3D=$4
INPUTBASE=$(basename $INPUT)
EXEC=generateTree2D
EXEC3D=generateTree3D

    




if test -f "$INPUT"
then
  convert ${INPUT} input.pgm
  ${EXEC} -n ${NBTERM} -a ${APERF}  -d input.pgm >> algo_info.txt
  type ${EXEC}
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
  # Strangely enough we have to change with execution mode...
  chmod +x /srv/home/ipol/ipolDevel/ipol_demo/modules/demorunner/binaries/5555531082024/bin/generateTree3D

  ${EXEC3D} -n ${NBTERM} -a ${APERF}  -o result.obj
  cat stderr.txt
  key=$(basename $(pwd))
  demo_id=$(basename $(dirname $(pwd)))
  viewer_url="https://3dviewer.net#https://ipolcore.ipol.im/api/core/shared_folder/run/${demo_id}/${key}/result.obj,https://ipolcore.ipol.im/api/core/shared_folder/run/${demo_id}/${key}/result.mtl"
  iframe="<iframe id='3dviewerplayer' type='text/html' width='620' height='460' src='$viewer_url' "
  iframe="$iframe frameborder='5' scrolling='no' allowfullscreen webkitallowfullscreen mozallowfullscreen></iframe>"
  echo "url=$iframe" >> algo_info.txt
  echo "algoDim=3" >> algo_info.txt 
fi




