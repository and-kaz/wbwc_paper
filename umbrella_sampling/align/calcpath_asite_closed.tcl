####                                                              ####
#### 	Calculates pathCV and hb variable for the frames 	  ####
#### 	outputted in the previous script, and writes it into      ####
####	a file for further pipeline (also writes atom indices)    ####
######################################################################
## 0. load frames and images (
set reacname "GU_wbwc_neb_full"
set sysname "asite_closed"
set iseq "neb_nodpt"

set psfloc    "../init/${sysname}.psf"
set framefold "../systems/frames/${sysname}/${reacname}"
set imagefold "../systems/images/${sysname}_${reacname}_${iseq}"
####
## Frames
mol new $psfloc
foreach f [lsort [glob $framefold/*.pdb]] {
mol addfile $f waitfor all
}
##
## Images
mol new $psfloc
foreach f [lsort [glob $imagefold/*.pdb]] {
mol addfile $f waitfor all
}
#################
## select path atoms in frames based on images
set pathselim [atomselect 1 "beta 1"]
set pathatlist [$pathselim get index]
set pathselframe [atomselect 0 "index $pathatlist"]

### Calculate path
#### SET LAMBDA MANUALLY !!!
####################
set lambdaval 300
#################
set colvardat [open "../colvars/cv_${sysname}_${reacname}.dat" w]

##  measure O6-H3 (assumes no H3 tautomer of Gx)
## since frames don't have BP beta'd, select atoms in images, then get indexes to select them in frames)

###!! if starting from G[enol], names of the atoms are different!!! ###
set o6im [atomselect 1 "(same residue as beta 1) and name O6"]
set h3im [atomselect 1 "(same residue as beta 1) and resname GUA and name H6"]
set n1im [atomselect 1 "(same residue as beta 1) and resname GUA and name N1"]
set h1im [atomselect 1 "(same residue as beta 1) and name H3"]

set o6ind [$o6im get index]
set h3ind [$h3im get index]
set n1ind [$n1im get index]
set h1ind [$h1im get index]

#puts "o6ind: ${o6ind}"

set o6 [atomselect 0 "index ${o6ind}"]
set h3 [atomselect 0 "index ${h3ind}"]
set n1 [atomselect 0 "index ${n1ind}"]
set h1 [atomselect 0 "index ${h1ind}"]

## print atomnumbers (index +1) of O6 and H3 for colvar definition
set o6n [expr ([$o6 get index] +1)]
set h3n [expr ([$h3 get index] +1)]
set n1n [expr ([$n1 get index] +1)]
set h1n [expr ([$h1 get index] +1)]
puts $colvardat "#O6: ${o6n}, H3: ${h3n}, N1: ${n1n}, H1: ${h1n}"

set nframes [molinfo 0 get numframes]
set nimages [molinfo 1 get numframes]

for {set i 0} {$i < $nframes} {incr i} {
  $pathselframe frame $i
  $o6 frame $i
  $h3 frame $i
  $n1 frame $i
  $h1 frame $i
  set framenumer 0
  set framedenom 0
  for {set j 0} {$j < $nimages} {incr j} {
    $pathselim frame $j
    set ijrmsd [measure rmsd $pathselframe $pathselim]
    set framenumer [expr {$framenumer + $j * exp(- $lambdaval * $ijrmsd * $ijrmsd)}]
    set framedenom [expr {$framedenom + exp(- $lambdaval * $ijrmsd * $ijrmsd)}]
  }
  set ipath [expr {1.0 / ($nimages - 1.0) * $framenumer / $framedenom}]
  set o6c [measure center $o6]
  set h3c [measure center $h3]
  set n1c [measure center $n1]
  set h1c [measure center $h1]
  set o6h3 [veclength [vecsub $h3c $o6c]]
  set n1h1 [veclength [vecsub $h1c $n1c]]
  set hb [expr {$o6h3 - $n1h1}]
  puts $colvardat "$ipath   $hb"
}

close $colvardat  
puts "done!"
    
