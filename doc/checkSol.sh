#!/bin/bash

export MYWDIR="/data0/devAlg/tests"

export DataList="collCheckList.txt"
export PREFDIR="./"
: ${ErrDet:=";GPU_rec_trd.minTrackPt=0.3;;GPU_rec_trd.trkltResRPhiIdeal=0.2;GPU_rec_trd.trkltResVsTanPhiMisalign=1.;AlignParams.extraClsErrYTRD=0.1;AlignParams.extraClsErrZTRD=1.;AlignParams.extraClsErrYTOF=0.2;AlignParams.extraClsErrZTOF=0.2;"}
for ii in $(seq 0 6) ; do
    ErrDet+="AlignParams.extraClsErrYITS[$ii]=0.0003;AlignParams.extraClsErrZITS[$ii]=0.0003;"
done
export EXTRACONF="AlignParams.debugTracksFraction=1;AlignParams.maxTracksPerMultiTrackPV=10;AlignParams.minPt=0.8;AlignParams.meanPtB0=0.7;AlignParams.corrType=1;AlignParams.maxChi2Ndf=14.f;$ErrDet;"
if [[ -d $(pwd)/GLO/Config/GeometryAligned ]] ; then
    export REMAP="--condition-remap file://$(pwd)=GLO/Config/GeometryAligned"
fi
export OutType="ProcessData"
export MAXTR=${MAXTF:-2000}

${MYWDIR}/badProcTrdTof.sh

