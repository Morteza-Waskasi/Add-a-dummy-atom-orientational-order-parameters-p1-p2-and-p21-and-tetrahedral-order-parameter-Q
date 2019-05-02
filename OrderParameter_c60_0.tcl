# This code calculates the orientational order parameters p1, p2 and p21 and tetrahydral order parameter Q.
# For more details see: 
# Sarhangi, Waskasi, Hashemianzadeh and Matyushov, Phys. Chem. Chem. Phys., (2018), DOI: 10.1039/c8cp05422c
# For simulations details see: 
# Sarhangi, Waskasi, Hashemianzadeh, Martin  and Matyushov,  J. Phys. Chem. C 2018, 122, 17080−17087

package require pbctools

if { $argc != 10 } {
   puts "This vmd script requires 9 arguments."
   puts "Please try again."
   exit
}
puts "Starting Analysis - date/time"
date

set state        [ lindex $argv 0 ] ;# The name of simulation, c60_0 in this case.
set dcdDir       [ lindex $argv 1 ] ;# Path to trajectories (dcd)
set psfFile      [ lindex $argv 2 ] ;# Path to psf file
set shellDst     [ lindex $argv 3 ] ;# The shell distance for first hydration shell
set ns           [ lindex $argv 4 ] ;# The number of nano second. 
set temp         [ lindex $argv 5 ] ;# Temperature 
set p1p2outDir   [ lindex $argv 6 ] ;# Output file for average values
set p1p2outDird  [ lindex $argv 7 ] ;# Output fiel for every single value (distributions)
set Qshell       [ lindex $argv 8 ]

set ofileNm ${p1p2outDir}/p1_p2_${state}_${temp}_${ns}ns.dat     ;# name of output file (average values)
set ofileNmd ${p1p2outDird}/p1_p2d_${state}_${temp}_${ns}ns.dat  ;# name of output file (distribution) 

set dcdFile ${state}NVT.dum.${ns}ns.dcd                          ;# name of dcd file


set traj [ mol load psf ${psfFile} dcd ${dcdDir}/${dcdFile}]      ;# loading psf and dcd 
pbc wrap -centersel "segname DU FU" -center com -compound res -all ;# move the C60(with dummy atom at the center) to the center of box
set nFrames [ molinfo $traj get numframes ]   ;# get total number of frames 
puts [ format "Reading %i frames." $nFrames ] ;# write total number of frames 
puts "Number of atoms in the box = [ molinfo $traj get numatoms ]" ;# write total number of atoms in the box
set ofile  [ open ${ofileNm} w ]   ;# open output file
set ofiled [ open ${ofileNmd} w ]  ;# open output file


for {set f 0} {$f < $nFrames } {incr f 1} {
   puts "********************"
   puts "** FRAME : $f     **"
   puts "********************"

   set wat1st [atomselect $traj " same residue as resname SPCE and within ${shellDst} of resname DUM" frame $f] ;# select waters within the shell distance 
   set watlstindex [ $wat1st get index ]               ;# get index of all atoms of water 
   set dum [atomselect $traj "resname DUM" frame $f]   ;# select dummy atom
   set rdum [ lindex [ $dum get {x y z} ] 0 ]          ;# get the coordinate of dummy atom
   puts -nonewline [ format "FRAME %i: " $f ]

   set numWatShellf [ expr ( [$wat1st num]/3 )  ]      ;# print the number of waters whithin the shell distance
   puts "Number of waters in first solvation shell for frame $f is [ expr ( [$wat1st num]/3 )  ]"
     
   set oxlst [ atomselect $traj "index $watlstindex and name OH2" frame $f ]  ;# selecting  oxygens in the first shell
   set oxlstindex [ $oxlst get index ]                                        ;# get indices of oxygens in the first shell
 
   set i 0
   set sump1fi 0
   set sump1sfi 0
   set sump2fi 0
   set sump2sfi 0
   set sump21fi 0
   set sump21sfi 0
   set sumHbondNumber 0
   set sumTdQ  0
   set p2 0

############### ******** H-Bonds *************  ########################
# calculates number of H-Bonds in the shell with distance ${shellDst}
set hbond [measure hbonds 3.5 30 $wat1st]
set   hbondNumber [llength [lindex $hbond 0]]
########################################################################

    
      foreach i $oxlstindex {                                 ;# loop over oxygens in the first shell
        set oxi [ atomselect $traj "index $i" frame $f ]      ;# selecting the oxygen
        set oxires [ $oxi get resid ]                         ;# get resid of oxygen     
        set roxi  [ lindex [ $oxi get {x y z} ] 0 ]           ;# get coordinate of oxygen 
        set hy1i [ atomselect $traj "resid $oxires and name H1" frame $f ] ; # selecting hydrogen number 1 connected to the oxygen 
        set hy2i [ atomselect $traj "resid $oxires and name H2" frame $f ] ; # selecting  hydrogen number 2 connected to the oxygen
        set rhy1i [ lindex [ $hy1i get {x y z} ] 0 ]            ;# get coordinate of H1
        set rhy2i [ lindex [ $hy2i get {x y z} ] 0 ]            ;# get coordinate of H1
        set vechy1i    [ vecnorm [ vecsub $rhy1i $roxi ] ]      ;# get normal vector which connects oxygen to H1 
        set vechy2i    [ vecnorm [ vecsub $rhy2i $roxi ] ]      ;# get normal vector which connects oxygen to H2
        set wnormveci   [vecnorm [ vecadd $vechy1i $vechy2i]]   ;# normal vector represents direction of dipole moment
        set watnormi    [vecnorm [ veccross $vechy1i $vechy2i]] ;# normal vector which is normal to the plane of water
        set Rodum   [ vecnorm [ vecsub $roxi  $rdum ] ]         ;# normal vector which connects dummy atom to the oxygen
#        set normk   [vecnorm [veccross $wnormveci $watnormi]]
#        set normk   [vecnorm [veccross $Rodum $wnormveci ]]
#        set kDwnormi [vecnorm [vecdot $normk $watnormi ]]
        set costhetai [ vecdot $wnormveci $Rodum ]              ;# calculate the angle between water's dipole moment and
                                                                # the normal vector which connects dummy atom to the oxygen
                                                                # or cos(theta) which is p1, first order parameter
                                                                # see eq 13 in PCCP : Phys. Chem. Chem. Phys., (2018), DOI: 10.1039/c8cp05422c
#        set thetai  [expr acos($costhetai)]
#        set sinthetai [ veclength [veccross $wnormveci $Rodum ]]
        set sump1fi [expr $sump1fi + $costhetai]                 ;# sum all cos(theta)         
        set sump1sfi [expr $sump1sfi + ($costhetai*$costhetai)]  ;# sum of square of p1
        set p2    [expr (3*($costhetai*$costhetai)-1)*0.5]       ;# second order orientational parameter
        set sump2fi [expr $sump2fi + $p2 ]                       ;# sum p2
 
        set sump2sfi [expr $sump2sfi + $p2*$p2 ]                 ;# sum of square of p2
        set normnm [vecnorm [veccross $wnormveci $Rodum ]  ]
        set chi  [ expr acos( [vecdot $normnm $watnormi ])]      ;# angle chi which is the angle between the planes of the water
                                                                 # molecule and that of the dipole moment and the radial direction. 
                                                                 # see figure 1 in PCCP paper 
#        set p21fi [ expr ($sinthetai*$sinthetai)*cos(2* $chi)]
        set p21fi [ expr (0.5*(1-($costhetai)*($costhetai))) * cos(2* $chi)]  ;# calculates p21 see 
        set sump21fi [expr $sump21fi + $p21fi ]
        set sump21sfi [expr $sump21sfi + $p21fi*$p21fi]           ;# sum of square of p21      


########## **********calculation of tetrahedral order parameter Q*********  #######
################### see eq in PCCP                        #########################
                                  set   f4ox    {};
                                  set   nlist   {};
         # to calculate the tetrahydral order parameter, 4 nearst water molecuel should be selected.       
         # $Qshell is a shell with radius of  $Qshell A from oxygen with index i

        set Qwat [atomselect $traj " same residue as resname SPCE and within $Qshell of index $i" frame $f] ;# select all waters within $Qshell distance
        set Qwatlstindex [ $Qwat get index ]      ;# get index of waters within the $Qshell 
        set Qoxlst [ atomselect $traj "index $Qwatlstindex and name OH2" frame $f ] ;# selecting oxygens
        set Qoxlstindex [ $Qoxlst get index ]                                       ;# get index of oxygen



                      foreach j $Qoxlstindex {    ;# loop over all oxygens in Qshell 


                              set Qoxi [ atomselect $traj "index $j" frame $f ] ;# selecting oxygen 
                              set rQoxi  [ lindex [ $Qoxi get {x y z} ] 0 ]     ;# get coordinate of oxygen
                              set dist [veclength [vecsub $rQoxi $roxi ] ]      ;# calculate distance between the target oxygen and oxygens close to it
                              lappend nlist "$j $dist"                          
   
                              set f4ox [lrange [lsort -increasing -index 1 $nlist] 0 4 ] ;# select 4 nearst oxygens
 
                            $Qoxi delete 
}

                              set oxi1 [ atomselect $traj "index [lindex $f4ox 1 0]" frame $f ] ;# select O1 
                              set oxi2 [ atomselect $traj "index [lindex $f4ox 2 0]" frame $f ] ;# select O2
                              set oxi3 [ atomselect $traj "index [lindex $f4ox 3 0]" frame $f ] ;# select O3
                              set oxi4 [ atomselect $traj "index [lindex $f4ox 4 0]" frame $f ] ;# select O4

                              set roxi1 [ lindex [ $oxi1 get {x y z} ] 0 ] ;# get coordinate of  O1 
                              set roxi2 [ lindex [ $oxi2 get {x y z} ] 0 ] ;# get coordinate of  O2
                              set roxi3 [ lindex [ $oxi3 get {x y z} ] 0 ] ;# get coordinate of  O3
                              set roxi4 [ lindex [ $oxi4 get {x y z} ] 0 ] ;# get coordinate of  O4

                              set vecox1    [ vecnorm [ vecsub $roxi1 $roxi ] ] ;# get normal vector which connects the targent oxygen and O1
                              set vecox2    [ vecnorm [ vecsub $roxi2 $roxi ] ] ;# get normal vector which connects the targent oxygen and O2
                              set vecox3    [ vecnorm [ vecsub $roxi3 $roxi ] ] ;# get normal vector which connects the targent oxygen and O3
                              set vecox4    [ vecnorm [ vecsub $roxi4 $roxi ] ] ;# get normal vector which connects the targent oxygen and O4

                              set cosij1 [vecdot $vecox1 $vecox2 ] ;# calculate cos(theta), which theta is between $vecox1 $vecox2
                              set cosij2 [vecdot $vecox1 $vecox3 ] ;# calculate cos(theta), which theta is between $vecox1 $vecox3
                              set cosij3 [vecdot $vecox1 $vecox4 ] ;# calculate cos(theta), which theta is between $vecox1 $vecox4
                              set cosij4 [vecdot $vecox2 $vecox3 ] ;# calculate cos(theta), which theta is between $vecox2 $vecox3
                              set cosij5 [vecdot $vecox2 $vecox4 ] ;# calculate cos(theta), which theta is between $vecox2 $vecox4
                              set cosij6 [vecdot $vecox3 $vecox4 ] ;# calculate cos(theta), which theta is between $vecox3 $vecox4
         
                               # TdQ is tetrahydral order parameter see equation 12 in PCCP paper
                              set TdQ [expr 1 - (3.0/8.0)*( \
                                 ($cosij1 + 1/3.0)*($cosij1 + 1/3.0) +  \
                                 ($cosij2 + 1/3.0)*($cosij2 + 1/3.0) +  \
                                 ($cosij3 + 1/3.0)*($cosij3 + 1/3.0) +  \
                                 ($cosij4 + 1/3.0)*($cosij4 + 1/3.0) +  \
                                 ($cosij5 + 1/3.0)*($cosij5 + 1/3.0) +  \
                                 ($cosij6 + 1/3.0)*($cosij6 + 1/3.0) ) ]

 #prints values in the output file
  puts $ofiled [format  " %.0f\t%.0f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t"  $ns $f  $costhetai $p2  $p21fi $chi  $TdQ] 
 
                                  array unset roxi1; 
                                  array unset roxi2; 
                                  array unset roxi3; 
                                  array unset roxi4; 
                                  array unset nlist; 
                                  array unset hbondNumber; 
                                  array unset f4ox; 

                                  $Qwat   delete ; unset Qwat
                                  $oxi1   delete ; unset oxi1
                                  $oxi2   delete ; unset oxi2
                                  $oxi3   delete ; unset oxi3
                                  $oxi4   delete ; unset oxi4
                                  $Qoxlst delete ; unset Qoxlst
                   
          set sumTdQ [expr $sumTdQ + $TdQ]

        $oxi  delete
        $hy1i  delete
        $hy2i   delete
      } 
  set avep21fi [expr ( $sump21fi/$numWatShellf)]
  set avep21sfi [expr ( $sump21sfi/$numWatShellf)]
  set avep1fi [expr $sump1fi / $numWatShellf]
  set avep1sfi [expr $sump1sfi / $numWatShellf]
  set avep2fi [expr $sump2fi /(1* $numWatShellf)]
  set avep2sfi [expr $sump2sfi /(1* $numWatShellf)]


  set aveQTd [expr $sumTdQ / $numWatShellf]
 #prints average values in the output file
  puts $ofile [format  " %.0f\t%.0f\t%.0f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f"  $ns $f  $numWatShellf $hbondNumber  $avep1fi $avep1sfi  $avep2fi $avep2sfi $avep21fi $avep21sfi $aveQTd] 

  puts "cleaning up"
     $oxlst  delete
     $wat1st delete
     $dum delete
  }
close $ofile
close $ofiled
mol delete all
exit

    

