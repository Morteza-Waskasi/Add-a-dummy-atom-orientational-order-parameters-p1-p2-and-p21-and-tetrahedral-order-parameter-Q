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
set shell        [ lindex $argv 7 ]

set outShD  [expr ($shellDst + $shell ) ]

set ofileNm ${p1p2outDir}/p1_p2_${state}_${temp}_${ns}ns.dat

set dcdFile ${state}NVT.dum.${ns}ns.dcd


set traj [ mol load psf ${psfFile} dcd ${dcdDir}/${dcdFile}]
pbc wrap -centersel "segname DU FU" -center com -compound res -all
set nFrames [ molinfo $traj get numframes ]
puts [ format "Reading %i frames." $nFrames ]
puts "Number of atoms in the box = [ molinfo $traj get numatoms ]"
set ofile [ open ${ofileNm} w ]
#set sump1fi 0
for {set f 0} {$f < $nFrames } {incr f 1 } {
   puts "********************"
   puts "** FRAME : $f     **"
   puts "********************"
#pbc wrap -centersel "protein" -center com -compound res -all
#   set wat1st [atomselect $traj "resname SPCE and within ${shellDst} of resname DUM" frame $f]
#   set wat1st [atomselect $traj " same residue as resname SPCE and within ${shellDst} of resname DUM" frame $f]
set sump1fi 0

set wat1st [atomselect $traj "(same residue as resname SPCE and within ${outShD} of resname DUM ) and (same residue as resname SPCE and  not within ${shellDst} of resname DUM)" frame $f] 
   set watlstindex [ $wat1st get index ]
   set dum [atomselect $traj "resname DUM" frame $f]
   set rdum [ lindex [ $dum get {x y z} ] 0 ]
   puts -nonewline [ format "FRAME %i: " $f ]

   set numWatShellf [ expr ( [$wat1st num]/3 )  ]
   puts "Number of waters in first solvation shell for frame $f is [ expr ( [$wat1st num]/3 )  ]"
      if { $numWatShellf==0 } {
      set sump1fi 0 
      puts "Number of waters in first solvation shell for frame  $f is zero"       
       } else {

  
# Standard simple if with else clause
#if {some condition} {
#    some conditionally executed script.
#} else {
#    some script to execute if the condition is not satisfied.
#}

     # get indices of oxygens in the first shell
   set oxlst [ atomselect $traj "index $watlstindex and name OH2" frame $f ]
   set oxlstindex [ $oxlst get index ]
 
   set i 0
   set sump1fi 0
   set sump1sfi 0
   set sump2fi 0
   set sump2sfi 0
   set sump21fi 0
   set sump21sfi 0
   set sumHbondNumber 0
   set sumTdQ  0

    # loop on oxygens
      foreach i $oxlstindex {
#      incr i
        set oxi [ atomselect $traj "index $i" frame $f ]
puts "oxi in i $i"
        set oxires [ $oxi get resid ]
        set roxi  [ lindex [ $oxi get {x y z} ] 0 ]
        set hy1i [ atomselect $traj "resid $oxires and name H1" frame $f ]
        set hy2i [ atomselect $traj "resid $oxires and name H2" frame $f ]
        set rhy1i [ lindex [ $hy1i get {x y z} ] 0 ]
        set rhy2i [ lindex [ $hy2i get {x y z} ] 0 ]
        set vechy1i    [ vecnorm [ vecsub $rhy1i $roxi ] ]
        set vechy2i    [ vecnorm [ vecsub $rhy2i $roxi ] ]
        set wnormveci   [vecnorm [ vecadd $vechy1i $vechy2i]]
        set watnormi    [vecnorm [ veccross $vechy1i $vechy2i]]
        set Rodum   [ vecnorm [ vecsub $roxi  $rdum ] ]
        set costhetai [ vecdot $wnormveci $Rodum ]
        set sump1fi [expr $sump1fi + $costhetai]


        $oxi  delete
        $hy1i  delete
        $hy2i   delete
      } 
     $oxlst  delete

      }
     $wat1st delete
     $dum delete
  puts $ofile [format  "%.6f"  $sump1fi] 
#  puts $ofile [format  " %.0f   %.6f   %.0f   %.6f   %.6f   %.6f   %.6f"  $f $sump1fi $numWatShellf $avep1fi $sump2fi $avep2fi $avep21fi] 
  puts "cleaning up"
#     $oxlst  delete
#     $wat1st delete
#     $dum delete
#     $watlstindex delete
#    unset  $rdum 

# animate delete all
  }

close $ofile
mol delete all
exit

    

