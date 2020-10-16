#!/bin/bash
#SBATCH -A hhcankaz
#SBATCH -p standard96
#SBATCH -N 1
#SBATCH --tasks-per-node 96
#SBATCH --mem 32G
#SBATCH -t 11:59:00
#SBATCH -J 6gsk_nen_1
##########
## 
module load slurm
###
# to prevent VMD multithreading
VMDFORCECPUCOUNT=1
##
ussys="6gsk_nen"
##
## starts from one kf, creates tmux and sends commands to them, then proceeds to next kf
cd kf_05/run
##
declare -A excluded
for i in {1,3,5,6,20,22,23,30,31,33}; do
excluded[$i]=1
done
### 1. concatenate trajectories (and rename colvars.state.old -> colvars.state (namd is stupid and does not allow to use *.old in config)
for i in {0..34}; do
  if ! [[ ${excluded[$i]} ]]; then
    trajname="${ussys}_us_${i}.colvars.traj"
    dcdname="${ussys}_us_${i}.dcd"
    if [ ! -d "${i}/traj1" ]; then
      mkdir "${i}/traj1"
      mkdir "${i}/traj2"
      mkdir "${i}/dcd"
    fi

    if [ -f "${i}/$trajname" ]; then
      mv ${i}/${ussys}_us_${i}.restart.colvars.state.old ${i}/${ussys}_us_${i}.restart.colvars.state
      if [ -f "${i}/traj1/${trajname}" ]; then
        cat "${i}/traj1/${trajname}" "${i}/$trajname" > "${i}/traj2/${trajname}"
        cp "${i}/traj2/${trajname}" "${i}/traj1/${trajname}"
        mv --backup=t "${i}/${dcdname}" "${i}/dcd/${dcdname}"
      else
        cp "${i}/$trajname" "${i}/traj1/"
        cp "${i}/$trajname" "${i}/traj2/"
        mv --backup=t "${i}/${dcdname}" "${i}/dcd/${dcdname}"
      fi
    fi
  fi
done
#######
## create tmux session, send run VMD (for restraints), and then NAMD commands to them
#######
for i in {0..34}; do
  if ! [[ ${excluded[$i]} ]]; then
    echo "starting kf05 ${i}"
    vmdres="vmdres_${i}.tcl"
    comm_vmd="/home/hhcankaz/bin/vmd -e ${i}/${vmdres}"
    comm_namd="namd2 ${i}/us_${i}.conf > ${i}/us_$i.log"
    tmux new-session -d -s "comm2_$i"
    tmux send-keys -t "comm2_$i" "$comm_vmd" C-m
    sleep 5
    tmux send-keys -t "comm2_$i" "$comm_namd" C-m
    sleep 1
  fi
done
#############################################
## next Kf #################################
############################################
sleep 6
######
cd ../../kf_005/run
##
declare -A excludedn
for i in {10..19}; do
excludedn[$i]=1
done
### 1. concatenate trajectories
for i in {0..34}; do
  if ! [[ ${excludedn[$i]} ]]; then
    trajname="${ussys}_us_${i}.colvars.traj"
    dcdname="${ussys}_us_${i}.dcd"
    if [ ! -d "${i}/traj1" ]; then
      mkdir "${i}/traj1"
      mkdir "${i}/traj2"
      mkdir "${i}/dcd"
    fi

    if [ -f "${i}/$trajname" ]; then
      mv ${i}/${ussys}_us_${i}.restart.colvars.state.old ${i}/${ussys}_us_${i}.restart.colvars.state
      if [ -f "${i}/traj1/${trajname}" ]; then
        cat "${i}/traj1/${trajname}" "${i}/$trajname" > "${i}/traj2/${trajname}"
        cp "${i}/traj2/${trajname}" "${i}/traj1/${trajname}"
        mv --backup=t "${i}/${dcdname}" "${i}/dcd/${dcdname}"
      else
        cp "${i}/$trajname" "${i}/traj1/"
        cp "${i}/$trajname" "${i}/traj2/"
        mv --backup=t "${i}/${dcdname}" "${i}/dcd/${dcdname}"
      fi
    fi
  fi
done
#######
## create tmux session, send run VMD (for restraints), and then NAMD commands to them
#######
for i in {0..34}; do
  if ! [[ ${excludedn[$i]} ]]; then
    echo "starting kf005 ${i}"
    vmdres="vmdres_${i}.tcl"
    comm_vmd="/home/hhcankaz/bin/vmd -e ${i}/${vmdres}"
    comm_namd="namd2 ${i}/us_${i}.conf > ${i}/us_$i.log"
    tmux new-session -d -s "comm3_$i"
    tmux send-keys -t "comm3_$i" "$comm_vmd" C-m
    sleep 5
    tmux send-keys -t "comm3_$i" "$comm_namd" C-m
    sleep 1
  fi
done
###########
## sleeping is required to keep the job running
###########
sleep 11h 55m





