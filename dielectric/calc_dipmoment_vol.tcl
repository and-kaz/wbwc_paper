##########
## script for VMD to parse MD trajectories to calculate volume and dipole moment for 
## the calculations of didlectric permittivity using the Kirkwood-Froehlich equations
#####
## to prepare a job input, use an example script to sed the UPPERCASE VARIABLES to the needed values
##
cd /home/kazantsev/MM/asite/STATE/AVERSION/reps/REPNUM
mol new /home/kazantsev/MM/asite/struct/STATE_AVERSION.psf
##
#load all DCDs from nvt2
foreach i [lsort -dictionary [glob -type d *]] {
if {[string is integer -strict $i]} {
mol addfile $i/STATE_AVERSION_REPNUM.dcd waitfor all
}
}
###############
## FUNCTIONS ##
###############
proc mean  L {
    expr ([join $L +])/[llength $L].
}
#
proc ladd l {
    if {![llength $l]} {return 0}
    return [expr [join $l +]]
}
#########
## 1. calc. the number of selections residues. choose the LARGEST SAMPLED SELECTION
## skip some eq.time (clumsy f-n, but saves memory by not storing a list)
proc get_asite_num {sel} {
  #set asite_num_dat [open "/home/kazantsev/MM/asite/parse/eps/asite_number_STATE_AVERSION_REPNUM_WATERCUTOFF_BPNUM.dat" w]
  set n [molinfo top get numframes]
  for {set i 0} {$i < $n} {incr i} {
    molinfo top set frame $i
    set asite [atomselect top "$sel"]
    set asite_ind [$asite get index]
    if {$i > 10} {
      if {[$asite num] > $asite_num} {
        set sel_ind [$asite get index]
        set sel_num [$asite num]
      }
    } else {
      set sel_ind [$asite get index]
      set sel_num [$asite num]
    }
    set asite_num [$asite num]
    $asite delete
  }
  puts "largest selection: $sel_num atoms"
  return $sel_ind
}
###########
## 2. calculate the number of molecules around selections, get the mean value
proc get_ref_wat {sel cut} {
  set n [molinfo top get numframes]
  set numlist []
  for {set i 0} {$i < $n} {incr i} {
    molinfo top set frame $i
    set wati [atomselect top "same residue as water and within $cut of $sel"]
    set numi [expr [$wati num]/3]
    lappend numlist $numi
  }
  ### stupid tcl doesn't even have a buil-in f-c for max/min... (mathfunc fails also)
  set numlistc [join $numlist ,]
  set nmax [expr max ($numlistc)]
  set nmin [expr min ($numlistc)]
  set nmean [expr round([mean $numlist])]
  puts "  --------- water molecules: ------------  "
  puts "min:${nmin},  max:${nmax},  mean:${nmean} "
  puts "  ---------------------------------------  "
  return $nmean
}
#####################
## 3. main f-n. for a reference number of wat molecules N, calculates N within cutoff of selection.
##    uses three (manually set) increments to switch between and avoid oscilations
####################
proc nwat {refnum sel cut volstep asite_ind} {
set cut_incr0 0.2
set cut_incr1 0.03
set cut_incr2_0 0.01
##################
molinfo top set frame 0
set asite_const [atomselect top "index $asite_ind"]
##
set asite_charge [expr round([ladd [$asite_const get charge]])]
###
set log_nwat [open "/home/kazantsev/MM/asite/scripts/epsilon/STATE_AVERSION/logs/STATE_AVERSION_REPNUM_WATERCUTOFF_BPNUM.log" w]
set out_dip_nwat [open "/home/kazantsev/MM/asite/parse/eps/dipole/asite_dip_STATE_AVERSION_REPNUM_WATERCUTOFF_BPNUM.dat" w]
file mkdir "/home/kazantsev/MM/asite/parse/eps/volume/STATE_AVERSION_REPNUM_WATERCUTOFF_BPNUM"
#######
set n [molinfo top get numframes]
puts $log_nwat "analyzing STATE_AVERSION, REPNUM, BPNUM bp with WATERCUTOFF cutoff..."
puts $log_nwat "traj length: $n"
puts $log_nwat " ------------------- "
set stepcount 0
for {set i 0} {$i < $n} {incr i} {
  set cut_incr2 $cut_incr2_0
  molinfo top set frame $i
  set cuti $cut
  set wati [atomselect top "same residue as water and within $cuti of $sel"]
  set wat_ind [$wati get index]
  set numi [expr [$wati num]/3]
  
  set cut_conv []
  set convlen [llength $cut_conv]
  set flag_osc 0
  while {$numi != $refnum} {
    set convlen [llength $cut_conv]
      if {($convlen > 200) || ($flag_osc == 1) } {
        set cut_incr [expr $cut_incr2/($convlen/2)]
      } else {
    
      if {abs([expr $refnum - $numi]) < 8} {
        set convlen [llength $cut_conv]
        if {abs([expr $refnum - $numi]) < 2} {
          if {$convlen > 2} {
            if {([lindex $cut_conv end]==[lindex $cut_conv end-2]) && ([lindex $cut_conv end]!=[lindex $cut_conv end-1])} {
              set cut_incr [expr $cut_incr2/($convlen/2)]
              set flag_osc 1
            } else {
              set cut_incr $cut_incr2
            }
          } else {
            set cut_incr $cut_incr2
          }
        } else {
          set convlen [llength $cut_conv]
          if {$convlen > 2} {
            if {([lindex $cut_conv end]==[lindex $cut_conv end-2]) && ([lindex $cut_conv end]!=[lindex $cut_conv end-1])} {
              set cut_incr [expr $cut_incr1/($convlen/2)]
              set flag_osc 1
            } else {
              set cut_incr $cut_incr1
            }
          } else {
            set cut_incr $cut_incr1
          }
        }
      } else {
         set convlen [llength $cut_conv]
         if {$convlen > 2} {
            if {([lindex $cut_conv end]==[lindex $cut_conv end-2]) && ([lindex $cut_conv end]!=[lindex $cut_conv end-1])} {
              set cut_incr [expr $cut_incr0/($convlen/2)]
            } else {
              set cut_incr $cut_incr0
            }
          } else {
            set cut_incr $cut_incr0
          }
        }
    }
    if {$numi < $refnum} {
      set cuti [expr $cuti + $cut_incr]
    } else {
      set cuti [expr $cuti - $cut_incr]
    }
    lappend cut_conv $cuti
    ## DEBUG::
    #set convlen [llength $cut_conv]
    #if {$convlen > 1} {
    #  puts $log_nwat "convlen: $convlen,  $numi != $refnum, cutoff: $cuti"
    #}
    ### :: debug
    ## as a last resort: select the needed amount (basically exclude one extra closest mol. manually)
    if {$convlen > 300} {
      set wati [atomselect top "same residue as water and within $cuti of $sel"]
      set numi [expr [$wati num]/3]
      if {$numi > $refnum} {
        set ref_ind [expr $refnum * 3]
        set wat_ind0 [$wati get index]
        set wat_ind [lrange $wat_ind0 0 $ref_ind]
        puts $log_nwat "manual water picking activated (at $numi != $refnum)"
        break
      } else {
        set wat_ind [$wati get index]
      }
    } else {
      set wati [atomselect top "same residue as water and within $cuti of $sel"]
      set wat_ind [$wati get index]
      set numi [expr [$wati num]/3]    
    }  
    $wati delete
  }
  puts $log_nwat "$i: water:: $convlen  "
  ## for log files:
  ##set wat_radius $cuti
  ## now everything is converged. select the whole a-site and measure stuff
  #####
  set all_asite [atomselect top "index $asite_ind $wat_ind"]
  ##
  set dipoles [measure dipole $all_asite -debye -masscenter]
  set dipole [veclength $dipoles]
  puts $out_dip_nwat "$dipole"

  incr stepcount
  if {$stepcount == $volstep} {
    volmap density $all_asite -res 1.0 -radscale 1.0 -weight mass -o "/home/kazantsev/MM/asite/parse/eps/volume/STATE_AVERSION_REPNUM_WATERCUTOFF_BPNUM/vol_dens_f${i}.dx"
    set stepcount 0
  }
 $all_asite delete
}
puts $log_nwat "---------------------done-----------------------------------"
close $log_nwat
close $out_dip_nwat
}
#################
#################
### RUN #########
## get selections
set asite_selection "ASITESELECTION"
set water_center "WATERSELECTION"

## 1. finds the largest selection of the trajectory. this selection will be constant during nwat function

set asite_index [get_asite_num $asite_selection]

## 2. get the mean waters number for given radius. will be kept constant in nwat

set ref_water_number [get_ref_wat $water_center WATERCUTOFF]

## 3. run the nwat function

nwat $ref_water_number $water_center WATERCUTOFF ASITEVOLUMESTEP $asite_index

exit



