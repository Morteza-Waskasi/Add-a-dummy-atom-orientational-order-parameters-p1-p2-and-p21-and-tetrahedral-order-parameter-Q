puts "Starting Analysis - date/time"
date
# This code adds one (or more) dummy atom to dcd. 
# For simulations details see:
# Sarhangi, Waskasi, Hashemianzadeh, Martin  and Matyushov,  J. Phys. Chem. C 2018, 122, 17080−17087
# OR 
# Sarhangi, Waskasi, Hashemianzadeh and Matyushov, Phys. Chem. Chem. Phys., (2018), DOI: 10.1039/c8cp05422c 

# To create psf file for the system, which contains dummy atom, 
# I would suggest to add one dummy atom to the pdb and use psfgen to create psf file. 
# However, the parameters (name, resname , segname , mass ,resid) should be provided for psfgen.

                              ; # Three parameters feeded by bash script 
set ns    [ lindex $argv 0 ]  ; # The number of nano second. 
set ps    [ lindex $argv 1 ]  ; # The number of pico second. 
set typ2  [ lindex $argv 2 ]  ; # The type of system, c60_0, C60 molecule with charge of zero.

set styp  NVT   ;# Type of simulation, NVT in this case.
set temp  300   ;# Temperature 

set baseDir      /scratch/02616/mmw/add_dummy/c60/${typ2}/${temp}    ; # Path to base directory, or PWD
set pdbPsfDir    ${baseDir}/molPdbPsf                                ; # Path to psf file
set psfFile      ${typ2}.psf                                         ; # psf file for C60 molecule with charge zero 
set psfFileDUM   ${typ2}dum.psf                                      ; # psf file for C60 molecule plus dummy atom

set dcdDir    /scratch/02616/mmw/c60/${typ2}/${styp}/${temp}/${ns}ns/${ps}ps    ;# Path to trajectories (dcd) 
set dcdFile   ${typ2}NVT.dcd                                                    ;# Name of dcd file 
set dumdcdDir ${baseDir}/traj_dum/${ns}ns                                       ;# Path to save new dcd(original dcd + dummy atom)                                     

package require pbctools
package require topotools

set traj       [ mol load psf ${pdbPsfDir}/${psfFile} dcd ${dcdDir}/${dcdFile} ] ;# loading psf and dcd
set trajdum    [ mol load psf ${pdbPsfDir}/${psfFileDUM}]                        ;# loading psf file for C60 molecule plus dummy atom

set nf [molinfo $traj get numframes]         ; # get the number of frames in the dcd
# loop over all frames
for {set frame 0} {$frame < $nf} {incr frame} {
 animate goto $frame    ; # go to the target frame
 set numdummy 1         ; # the number of dummy atoms to be add to the original dcd 
 set dum [mol new atoms $numdummy] ;# mol new is used to create a new atom (in this case dummy atom)
 unset numdummy                    ;# unset to free memory 
 animate dup $dum                  ;# Duplicate the given frame (default ``now'') of molecule molId and add the new frame to this molecule.
 set dumy [atomselect $dum all]    ;# select added dummy atom
 $dumy set name DUM ; $dumy set resname DUM ;$dumy set segname DU;$dumy set mass 0.001 ; $dumy set resid 0 ;# giving physical properties to the atom
 $dumy set {x y z} {{1.30 1.20 1.10}}          ;# giving coordinate 
 unset dumy                                    ;# unset to free memory
 set mol [::TopoTools::mergemols "$traj $dum"] ;# mergemols combines multiple separate molecules into one file. This is non-destructive and will create a new molecule.
 set carbon [atomselect $mol "resname FUL"]    ;# Select the C60 molecule 
 set mergdumy [atomselect $mol "resname DUM"]  ;# select the created dummy atom
 $mergdumy set {x y z} {{0.30 0.20 0.10}}      ;# set coordinate for dummy atom 
 $mergdumy moveto [measure center $carbon weight mass] ;# move the dummy atom to the center of C60
 
 unset mergdumy  ;# unset to free memory
 unset carbon    ;# unset to free memory
 set sel [atomselect $mol all]  ;# select atoms in whole box
#animate write pdb $frame.pdb sel $sel $mol
 animate write dcd ${dumdcdDir}/${typ2}NVT.${ns}ns.$frame.dcd sel $sel $mol  ;# write every frame on hard drive
#animate write psf $frame.psf sel $sel $mol
 unset sel      ;# unset to free memory
 mol delete $dum    
 unset dum       ;# unset to free memory
 mol delete $mol
 unset mol     ;# unset to free memory 
 }
 

for {set frame 0} {$frame < $nf} {incr frame} {
animate read dcd ${dumdcdDir}/${typ2}NVT.${ns}ns.$frame.dcd waitfor all  ;# loading every singel dcd in to memeory 
}
pbc wrap -centersel "segname FU" -center com -compound res -all          ;# move the C60(with the dummy atom at the center) to the center of box
animate write dcd ${dumdcdDir}/${typ2}NVT.dum.${ns}ns.dcd waitfor all    ;# write all of the frames to one dcd

exit


