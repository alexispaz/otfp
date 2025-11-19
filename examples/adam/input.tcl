# extraer_coors.tcl
mol new ../systems/alad.psf
mol addfile ../OTFP2D_periodic/job0.dcd step 200 waitfor all 

set all [atomselect 0 "all"]
set n [molinfo top get numframes]
for {set i 0} {$i < $n} {incr i} {
    $all frame $i
    $all writepdb [format "input/frame_%05d.pdb" $i]
}
quit

