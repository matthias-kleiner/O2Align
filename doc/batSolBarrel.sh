#!/bin/bash

: ${MYWDIR:="/data0/devAlg/tests"}

runProc() # args: destDir var=value [var=value ...]
{
  if [[ $# -lt 2 ]]; then
    echo "runProc expects arguments: <destDir> <var=value> [<var=value> ...]"
    exit 1
  fi
  local destDir=$1
  shift
  for assignment in "$@"; do
    if [[ ! $assignment =~ ^[a-zA-Z_][a-zA-Z0-9_]*=.*$ ]]; then
      echo "runProc expects variable assignments in the form <var>=<value>: $assignment"
      exit 1
    fi
  done
  echo "Running alignment with assignments $* in dir. $destDir"
  [[ ! -d $destDir ]] && mkdir $destDir
  pushd $destDir
  env "$@" $MYWDIR/sol.sh
  popd
}

export STARTMILLE=${STARTMILLE:-"${MYWDIR}/best/IBEdge/"}
export DetList=${DetList:-"ITS-TRD,ITS-TOF,ITS-TRD-TOF"}

#export EXTRACONF="AlignParams.MVTimeSlotsJson=${MYWDIR}/../data/MV_B0BPBM.json;AlignParams.statInFile=${MYWDIR}/Mille_IBEdgeS4R/procStartColl.root;"
export EXTRACONF=";AlignParams.statInFile=${MYWDIR}/Mille_IBEdgeS4R/procStartColl.root;AlignParams.corrType=1;"

#if [[ -n ${STARTMILLE} ]] &&  [[ ${STARTMILLE} != "none" ]] ; then
#  EXTRACONF+=";AlignParams.startFromMilleRes=${STARTMILLE}/millepede_labeled.res;"
#fi

export SKIPCHECKRESID=1


runProc ${MYWDIR}/solB0A1 Config=/home/shahoian/alice/O2Align/doc/ITSEnvTRDTOF_RB.json MPSTEER=${MYWDIR}/mpSteerBarrelV0.txt
#runProc ${MYWDIR}/solB0A1 Config=/home/shahoian/alice/O2Align/doc/ITSEnvTRDTOF_RB.json MPSTEER=${MYWDIR}/mpSteerBarrelV0TT0.txt

