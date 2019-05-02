package require pbctools
package require pbctools

if { $argc != 9 } {
   puts "This vmd script requires five arguments."
   puts "Please try again."
   exit
}
puts "Starting Analysis - date/time"
date

set state        [ lindex $argv 0 ]
set dcdDir       [ lindex $argv 1 ]
set psfFile      [ lindex $argv 2 ]
set shellDst     [ lindex $argv 3 ]
set ns           [ lindex $argv 4 ]
set temp         [ lindex $argv 5 ]
set p1p2outDir   [ lindex $argv 6 ]
#set p1p2outDird  [ lindex $argv 7 ]
set Qshell       [ lindex $argv 7 ]

set ofileNm ${p1p2outDir}/p1_p2_${state}_${temp}_${ns}ns.dat
#set ofileNmd ${p1p2outDird}/p1_p2d_${state}_${temp}_${ns}ns.dat

set dcdFile ${state}NVT.dum.${ns}ns.dcd


set traj [ mol load psf ${psfFile} dcd ${dcdDir}/${dcdFile}]
pbc wrap -centersel "segname DU FU" -center com -compound res -all
set nFrames [ molinfo $traj get numframes ]
puts [ format "Reading %i frames." $nFrames ]
puts "Number of atoms in the box = [ molinfo $traj get numatoms ]"
set ofile  [ open ${ofileNm} w ]
#set ofiled [ open ${ofileNmd} w ]

for {set f 0} {$f < $nFrames } {incr f 10} {
   puts "********************"
   puts "** FRAME : $f     **"
   puts "********************"
#pbc wrap -centersel "protein" -center com -compound res -all
#   set wat1st [atomselect $traj "resname SPCE and within ${shellDst} of resname DUM" frame $f]
   set wat1st  [atomselect $traj " same residue as resname SPCE and within ${shellDst} of resname DUM" frame $f]
   set wat2st  [atomselect $traj " same residue as resname SPCE and not within ${shellDst} of resname DUM" frame $f]
   set watDist [atomselect $traj " same residue as resname SPCE and within ${Qshell} of resname DUM" frame $f]
   set watall  [atomselect $traj " same residue as resname SPCE " frame $f]
   puts -nonewline [ format "FRAME %i: " $f ]

   set numWatShellf [ expr ( [$wat1st num]/3 )  ]
   set numWatShellfDist [ expr ( [$watDist num]/3 )  ]
   puts "Number of waters in first solvation shell for frame $f is [ expr ( [$wat1st num]/3 )  ]"
 
   set HbondNumber1   0
   set HbondNumber12  0
   set HbondNumberAll 0
   set HB12 0
############### ******** H-Bonds *************  ########################
set hbond [measure hbonds 3.5 30 $wat1st]
set   hbondNumber1 [llength [lindex $hbond 0]]
########################################################################
set hbondDist [measure hbonds 3.5 30 $watDist]
set   hbondNumberDist [llength [lindex $hbondDist 0]]
###################################################################
set hbond12 [measure hbonds 3.5 30 $wat1st $wat2st]
set   hbondNumber12 [llength [lindex $hbond12 0]]
########################################################################
set hbondAll [measure hbonds 3.5 30 $watall]
set   hbondNumberAll [llength [lindex $hbondAll 0]]

set HB12 [expr ($hbondNumber1 + $hbondNumber12) ]

#  puts $ofiled [format  " %.0f\t%.0f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t"  $ns $f  $costhetai $p2  $p21fi $chi  $TdQ] 
 


  puts $ofile [format  " %.0f\t%.0f\t%.0f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t"  $ns $f  $numWatShellf $hbondNumber1 $hbondNumber12  $HB12 $numWatShellfDist $hbondNumberDist $hbondNumberAll ] 
#  puts $ofile [format  " %.0f   %.6f   %.0f   %.6f   %.6f   %.6f   %.6f"  $f $sump1fi $numWatShellf $avep1fi $sump2fi $avep2fi $avep21fi] 
  puts "cleaning up"
  }

close $ofile
mol delete all
exit

    

