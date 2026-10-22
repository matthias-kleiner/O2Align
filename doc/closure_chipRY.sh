#!/bin/bash
# Closure test of the IB chip RY (rotation about the sensor normal): the sawtooth of dy vs Z with the 3 cm
# period seen in solA2v5 has a slope equal to the fitted chip RY, which asks whether the chain
# derivative -> pede -> AlignParam -> geometry is consistent for RY.
#
# The built-in misAlgJson injection cannot test it: it is ITS3-only (6 sensors) and it reuses the fit's own
# derivative chain, so a sign/pivot error in that chain would cancel. Here the misalignment is injected in the
# GEOMETRY (macro/injectChipRY.C), the data are untouched.
#
# Cases, both started from the geometry of the solution under test (BASE/ITS_alignment.root):
#   ctl : no injection. Self-consistent chain => the fit returns chip RY ~ 0 (the sawtooth must not come back).
#         A fit returning ~ -(its previous RY) means a sign error, ~ same RY again means RY is not applied.
#   inj : IB chip RY injected (+600, -600, +300 urad per layer, + 150 urad gaussian per chip).
#         Consistent chain => fitted RY = -injected (slope -1 in chipRYSummary.py).
#
# Usage (nothing is run until you call it; stages 1+3 are the usual procAll.sh ones):
#   ./closure_chipRY.sh ctl | inj | summary
#   WORK=/data0/devAlg/tests/closureRY BASE=/data0/devAlg/tests/solA2v5 ./closure_chipRY.sh inj
# Run it in the O2 + Millepede-II environment (alienv enter O2PDPSuite/latest Millepede-II/latest).
#
# ASSUMPTIONS to check against how solA2v5 was produced:
#  - ConfString below reproduces its options (dofConfigJson = ITS_RB_v5.json, procAll.sh defaults otherwise);
#  - one TF (the procAll.sh default DataDir) is enough for the layer-mean RY (12*9...20*9 chips per layer);
#    set MAXTF / DataDir for more;
#  - the initial alignment is picked from <case>/ITS/Calib/Align/snapshot.root by --condition-remap:
#    check "Building MP2AlignParams for detector ITS on top of initial alignment" in the stage 3 log.

set -e
CASE=${1:?usage: $0 ctl|inj|summary}
: ${BASE:=/data0/devAlg/tests/solA2v5}
: ${WORK:=/data0/devAlg/tests/closureRY}
: ${REPO:=/home/shahoian/alice/O2Align}
: ${Config:=$REPO/doc/ITS_RB_v5.json}   # ITS_RB_v6.json for the IB stave-RY variant
: ${ConfString:="AlignParams.dofConfigJson=$Config;AlignParams.maxChi2Ndf=15;AlignParams.useMultiTrackPVConstraint=2;AlignParams.maxTracksPerMultiTrackPV=20;AlignParams.minITSClsTrack=6;AlignParams.MVTimeSlotsJson=/data0/devAlg/data/MV_B0BPBM.json;"}
: ${DataDir:="/data0/devAlg/data/B0"}
PROC=/data0/devAlg/tests/proc0.sh

if [[ $CASE == summary ]]; then
  for c in ctl inj; do
    [[ -f $WORK/$c/result.json ]] || continue
    echo "== $c"
    python3 $REPO/macro/chipRYSummary.py $WORK/$c/result.json $( [[ $c == inj ]] && echo $WORK/inj/injectChipRY.txt )
  done
  exit 0
fi

D=$WORK/$CASE
mkdir -p $D/ITS/Calib/Align
cd $D

# 0. initial geometry of the case
if [[ $CASE == inj ]]; then
  root -l -b -q "$REPO/macro/injectChipRY.C(\"$BASE/ITS_alignment.root\", \"$D\", 6e-4, -6e-4, 3e-4, 1.5e-4)"
  mv ITS/Calib/Align/snapshot.root ./ITS_alignment.root
else
  cp $BASE/ITS_alignment.root ./     #$D/ITS/Calib/Align/snapshot.root
fi

TIMESTAMP="1746914681240"  ~/alice/mymacrosO2/makeCDB.sh

export Config ConfString
REMAP="--condition-remap file://$D=ITS/Calib/Align,GLO/Config/GeometryAligned"

# 1. Millepede data + steering
OutType="MilleData,MilleSteer" DataDir=$DataDir MAXTF=20 runOPTL="-b --run $REMAP" Lanes=${Lanes:-1} Nthreads=${Nthreads:-5} $PROC | tee stage1.log

# 2. pede, on the data of this case only
cp ../mpSteer.txt ./
###sed -n '1,/CFiles/p' $BASE/mpSteer.txt > mpSteer.txt
ls $D/mp2data.bin >> mpSteer.txt
pede mpSteer.txt | tee pede.log

# 3. result.json (+ ITS_alignment.root) from millepede.res
OutType="MilleRes" DataDir=$DataDir runOPTL="-b --run $REMAP" $PROC | tee stage3.log

python3 $REPO/macro/chipRYSummary.py result.json $( [[ $CASE == inj ]] && echo injectChipRY.txt )
