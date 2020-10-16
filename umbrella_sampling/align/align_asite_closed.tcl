####                                                              ####
####  I.  Fits a pdb structure to the reference trajectory        ####
####      in Gx-Ux base pair with n frames;                       ####
####  II. Outputs n frames and m<n images for pathCV              ####
####      using a selected predefined scheme, + rmsd matrix       ####                                                                        ####
####  
### Assumes the planar orientation of the base pair       ####
### otherwise alignment works bad and bonds are distorted ####
###
### Assumes only U, T and mnm5s2U (and enols) as Ux    ####
###
### When run, will save to disk (m+n) * size(sys.pdb) (could be a lot!)
#######################################################################
## 0. define the reference reaction trajectory and system to align
set reac "../ref/GU_wbwc_neb_full.pdb"
set sys "../init/asite_closed.pdb"
## select the predefined sequence of images for pathCV
set iseq "neb_nodpt"

set systype "asite"
##
set sysname [file rootname [file tail $sys]]
set reacname [file rootname [file tail $reac]]
###
## 1. loads the trajectory of the reference reaction, 
## finds its n frames, than loads the target pdb n times

mol new $reac waitfor all
set n [molinfo 0 get numframes]

mol new $sys
for {set i 1} {$i < $n} {incr i} {
mol addfile $sys
}

###
## 1. Initial alignement (of xyz ref to the pdb) ##
###
set refall [atomselect 0 all]
set qmall [atomselect 1 "beta 1"]

##Create arrays for indices
set qmset {}
set refset {}

### first alignent only by Gua
### aligning by both residues would create fluctuations in coordinates

foreach i [$qmall get index] {
  set iatom [atomselect 1 "index ${i}"]
  foreach j [$refall get index] {
    set jatom [atomselect 0 "index ${j}"]
    if {(([$iatom get resname] == [$jatom get resname]) &&
        ([$iatom get name] == [$jatom get name]) ||
	(([$iatom get name] == "H4") && ([$jatom get name] == "H3"))) ||
        (((([$iatom get resname] == "URA") && ([$jatom get resname] == "U8U")) || (([$jatom get resname] == "URA") && ([$iatom get resname] =="U8U"))) && 
        (([$iatom get name] == [$jatom get name]) ||
        (([$iatom get name] == "S2") && ([$jatom get name] == "O2")))) ||
        (((([$iatom get resname] == "URA") && ([$jatom get resname] == "THY")) || (([$jatom get resname] == "URA") && ([$iatom get resname] =="THY"))) &&
        (([$iatom get name] == [$jatom get name])))} {
        lappend qmset $i
        lappend refset $j
    }
    $jatom delete
  }
  $iatom delete
}

## Check of proper atom list
if {[llength $qmset] == [llength $refset]} {
set setlen [llength $refset]
} else {
puts "problem with common atom set finding!"
puts "target set length: [llength ${qmset}]"
puts "reference set length: [llength ${refset}]"
puts "check atom naming"
return -level 0 -code break
}

set qm [atomselect 1 "index ${qmset}"]
set ref [atomselect 0 "index ${refset}"]

##Creating order of indices in refset

set refsort [lsort -dictionary $refset]
set reforder {}

foreach i $refset {
  set o [lsearch $refsort $i]
  lappend reforder $o
}

##Fitting

for {set i 0} {$i < $n} {incr i} {
   molinfo 1 set frame $i
   molinfo 0 set frame $i
   set M0 [measure fit $ref $qm order $reforder] 
   $refall move $M0
}

## alignment of the reference reaction (as a rigid body) to the initial coord is done. now align initial coord to the reaction
#######################
## II. Atomwise alignment of the target to the trajectory ##

for {set j 0} {$j < $n} {incr j} {
   molinfo 1 set frame $j
   molinfo 0 set frame $j
   for {set i 0} {$i < $setlen} {incr i} {
     set targatom [atomselect 1 "index [lindex $qmset $i]"]
     set refatom [atomselect 0 "index [lindex $refset $i]"]
     set M1 [measure fit $targatom $refatom]
     if {[$targatom get name] == "C5"} {
       if {[$targatom get resname] == "URA"} {
         set 5group [atomselect 1 "resname URA and beta 1 and name C5 H5"]
         $5group move $M1
       } elseif {[$targatom get resname] == "THY"} {
         set 5group [atomselect 1 "resname THY and beta 1 and name C5 C5M H51 H52 H53"]
         $5group move $M1
       } elseif {[$targatom get resname] == "U8U"} {
           set 5group [atomselect 1 "resname U8U and beta 1 and name C5 C HC1 HC2 N HN1 HN2 CA HA1 HA2 HA3"]
           $5group move $M1
       } else {
           $targatom move $M1
         }
     } else {
         $targatom move $M1
     }
  }
}

puts "done with alignment"

## Move Ux's sugar along with the base (sugar is not in the QM, so it can be minimized easily)
## assumed no GUA modifications
set urabase [atomselect top "(beta 1) and not resname GUA"]
set refbasepos [measure center $urabase weight mass]
set urasugar [atomselect top "(same residue as (beta 1 and not resname GUA)) and not beta 1"]

for {set i 0} {$i < $n} {incr i} {
  molinfo top set frame $i
  set posnew [measure center $urabase weight mass]
  set vecmov [vecsub $posnew $refbasepos]
  $urasugar moveby $vecmov
}

### III. output frames and images
##
## 1. FOR ASITE MODELS, frames must contain restraints, which also use beta labels.
#### so beta label is added to the outer shell, but then again deleted before saving images
#### beta from the base pair must be removed for this, but then restored (stupid workaround)

file mkdir "../systems/frames/${sysname}/${reacname}"
file mkdir "../systems/images/${sysname}_${reacname}_${iseq}"
set rmsd_matrix [open "../systems/rmsmatrix/rmsd_matrix_${sysname}_${reacname}_${iseq}.dat" w]

## writing frames

set all [atomselect top all]
## ONLY WITH ASITE:
if {$systype == "asite"} {
  set bp_beta [atomselect top "beta 1"]
  set shell [atomselect top "same residue as (not (water or ions) and not within 27 of (chain V and resid 35))"]
  $all set beta 0
  $shell set beta 3
}

set prefix "../systems/frames/${sysname}/${reacname}/us_"

for {set i 0} {$i < $n} {incr i} {
molinfo top set frame $i
$all writepdb $prefix[format "%02i" $i].pdb
}

## ONLY WITH ASITE:
if {$systype == "asite"} {
  $all set beta 0
  $bp_beta set beta 1
}
######
puts "done with frames"
#####
## II.3 delete some frames (and beta labeling) to make a sequence of images for pathCV
## ALWAYS removes the Umod C5 group from the set
## workaround to remove it, but keep C5 (assumes no Gx!!!!):
set c5a [atomselect top "beta 1 and (not resname GUA) and name C5"]

if {$iseq == "neb_nodpt"} {
	$5group set beta 0
	$c5a set beta 1
	set hydrogens [atomselect top "beta 1 and name H3 H1 H21 H22 H6 H8"]
	$hydrogens set beta 0
	animate delete beg 25 end -1 
} else {
puts "select from the list of predefined sequences"
return -level 0 -code break
}

set n [molinfo top get numframes]
set bp0 [atomselect top "beta 1"]
set bp1 [atomselect top "beta 1"]

## write rmsd marix for analysis and lambda calculation

for {set i 0} {$i < $n} {incr i} {
  $bp0 frame $i
  set irmsd {}
  for {set j 0} {$j < $n} {incr j} {
    $bp1 frame $j
    set ijrmsd [measure rmsd $bp0 $bp1]
    lappend irmsd $ijrmsd
  }
  puts $rmsd_matrix $irmsd
}

close $rmsd_matrix


### writing the images

set all [atomselect top all]

set prefix "../systems/images/${sysname}_${reacname}_${iseq}/img_"

for {set i 0} {$i < $n} {incr i} {
molinfo top set frame $i
$all writepdb $prefix[format "%02i" $i].pdb
}

puts "done with images"


