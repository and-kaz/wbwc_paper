#!/bin/bash
#####
## example job input to run calc_dipmoment_vol.tcl on the trajectories of closed asite model (6gsk_nen)
## in this example, 3 MD replicas ca. 100 ns each are stored in parts of ca. 6 ns in separate folders
#####
state="6gsk"
aversion="nen"
watercut0=12
volume_step=2000
###
declare -A bpy
bpy=(["1bp"]=19 ["2bp"]=20 ["3bp"]=21)
declare -A bpv
bpv=(["1bp"]=36 ["2bp"]=35 ["3bp"]=34)
#######
#######
cd /home/kazantsev/MM/asite/scripts/epsilon/${state}_${aversion}/
######
for b in 1bp 2bp 3bp; do
  for i in rep1 rep2 rep3; do
    sely="${bpy[$b]}"
    selv="${bpv[$b]}"
    wat_center="chain Y and resid ${sely} and name C2"
    asite_sel="same residue as ((protein or nucleic) and not ((chain V and resid ${selv}) or (chain Y and resid ${sely}))) and within ${watercut0} of (chain Y and resid ${sely} and name C2)"
    sed -e "s/STATE/${state}/g" -e "s/AVERSION/${aversion}/g" -e "s/REPNUM/${i}/g" -e "s/BPNUM/${b}/g" -e "s/ASITESELECTION/${asite_sel}/g" -e "s/WATERSELECTION/${wat_center}/g" -e "s/WATERCUTOFF/${watercut0}/g" -e "s/ASITEVOLUMESTEP/${volume_step}/g" dipvol_asite_reps_ws_noions.tcl > dipvol_${state}_${aversion}_${i}_${b}_${watercut0}.tcl
    /home/kazantsev/app/vmd/bin/vmd -dispdev text -e dipvol_${state}_${aversion}_${i}_${b}_${watercut0}.tcl
    sleep 3
  done
done

 

