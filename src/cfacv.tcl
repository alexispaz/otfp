# tcl library for CFACV
# Cameron F Abrams 2009-14
# Sergio A Paz 2014-15

# This script is intended to be sourced before namd tclforce feature.  Therefore
# some of the procedures or variables used here are intrinsec of NAMD or can be
# found in the NAMD config file or in NAMD tcl master script for replica
# exechange.

# load the C-module
if {![info exists CFACV_BASEDIR]} {set CFACV_BASEDIR $::env(CFACV_BASEDIR)}
load ${CFACV_BASEDIR}/lib/libcfacv.so cfa_cvlibc


# load tcl procedures
source ${CFACV_BASEDIR}/include/data.tcl
source ${CFACV_BASEDIR}/include/restraints.tcl
source ${CFACV_BASEDIR}/include/cvs.tcl
source ${CFACV_BASEDIR}/include/chapeau.tcl

set SIG(TERM) 0

# BEGIN FUNCTION DEFINITION ###############################

proc cfacv_banner { argv } {
    global CFACV_VERSION
    global CFACV_BASEDIR
    print "================================================"
    print "CFACV) Collective Variables in VMD/NAMD v. $CFACV_VERSION"
    cfacvBanner
    print "CFACV) 2009-2014, Cameron F Abrams, Drexel University"
    print "CFACV) 2014-2015, Sergio A Paz, Drexel University"
    print "CFACV) 2016-2019, Sergio A Paz, Universidad Nacional de Córdoba - INFIQC/CONICET"
    print "================================================"
    print "Base directory: $CFACV_BASEDIR"
    print "argv: $argv"
    print "================================================"
    flush stdout
}

proc getArg { argv key altkey def } {
    set i 0
    while {$i < [llength $argv] && [lindex $argv $i] != "-$key" && [lindex $argv $i] != "-$altkey"} {
#	print "DB: getArg considers [lindex $argv $i]"
	incr i
    }
    if {$i < [llength $argv] && ([lindex $argv $i] == "-$key" || [lindex $argv $i] == "-$altkey")} {
	incr i
	set rv {}
	for {} {$i < [llength $argv] && [string index [lindex $argv $i] 0] != "-"} {incr i} {
	    lappend rv [lindex $argv $i]
#	    print "DB: getArg building rv as $rv"
	}
	return $rv
    }
    return $def
} 


proc Tcl_UpdateDataSpace { ds lC groups first timestep bias } {
    upvar $lC p
    global SIG

    # Move group center positions to dataspace
    set i 0
    foreach g $groups {
      ListToArray_Data [DataSpace_centerPosition $ds $i] $p($g)
      incr i
    }

    # Compute CV's within dataspace
    set localbias [DataSpace_ComputeCVs $ds]
    set bias [expr $localbias+$bias]
    # At this point, the (x,y,z) position data is no longer needed.
    # We can now write into its space the (x,y,z) restraint forces.

    # Compute restraining forces
    set rcode [DataSpace_RestrainingForces $ds $first $timestep $bias]
    if {$rcode == 1} {
       abort "Terminating due to internal OTFP signal."
    } elseif {$rcode == 2} {
       # abort "Terminating due to internal OTFP signal"
       if {$SIG(TERM) == 0 } {
         print "All restraints finish. Sending TERM signal."
         set SIG(TERM) 1
         # numsteps 10000
       }
    }

    # Transmit forces back to groups
    set i 0
    foreach g $groups {
        addforce $g [ArrayToList [DataSpace_centerForce $ds $i] 3]
        incr i
    }
    # puts "FUERZAS g1 [ArrayToList [DataSpace_centerForce $ds 0] 3]"

    # Add restraint energy to NAMD energy structure
    addenergy [DataSpace_RestraintEnergy $ds]
}

proc Tcl_ObserveDataSpace { ds cntrSel frame } {
    set i 0
    foreach sel $cntrSel {
	$sel frame $frame
#	print "DB: Tcl_ObserveDataSpace center $i [measure center $sel weight mass]"
	ListToArray_Data [DataSpace_centerPosition $ds $i] [measure center $sel weight mass]
	incr i
    }
    DataSpace_ComputeCVs $ds
}


# END FUNCTION DEFINITION ####################

# WARNING, this lines without a procedure
# implies that all variables used here are global

# Global variables
set CFACV_VERSION 0.30
set kB_kcm 0.001987191
set first 1
set masses {}

# Banner
cfacv_banner NAMD

# Set the pseudorandom number generator seed
if {![info exists seed]} {
  print "CFACV) setting seed from clock"
  set seed [clock clicks]
}
print "CFACV) random seed $seed"

# Read the template PDB file that identifies subdomain memberships
parse_pdb $labelPDB pdb
print "CFACV) nCenters $pdb(nctrs)"

# Allocate dataspace in C code
set ds [NewDataSpace $pdb(nctrs) $cv(num) $restr(num) $seed]
 
# Add CVs to the datapsace structure
cvs_setup pdb
 
# Add restraints to the dataspace structure
rlist_setup $cv(num)

# Add chapeau to the dataspace structure 
# TODO: improve chapeau global variable
if [info exists FES] {
  
    # Needed
    if {![info exists FES(max)]}    {error "Need FES(max)"}
    if {![info exists FES(min)]}    {error "Need FES(min)"} 
    if {![info exists FES(knots)]}  {error "Need FES(knots)"} 
    if {$restr(num) == 0} {error "No restraints to set a dataspace chapeau"}
   
    # Defaults
    if {![info exists FES(dm)]} {set FES(dm) 1}
    if {![info exists FES(period)]} {
      for {set i 0} {$i < $FES(dm)} {incr i} {lappend FES(period) 0}
    }

    # Replica exechange mode: 
    if {![info exists NUMREP]} {set NUMREP 1}
    # TODO: NUMREP podria ser un numero de chapeaus distintos (estilo el
    # problema sodio cloro) a obtener en una simulacion. En ese caso
    # corregir.

    # Saving for allocate chapeau functions
    chapeau_setup $NUMREP $ds $FES(dm) $FES(min) $FES(knots) $FES(max) $FES(period) $FES(start) $FES(solve.start) $FES(outfile) $FES(outfreq) $FES(outlvl) $FES(solve.freq)
}

# Some default values regarding output
if {![info exists TAMDoutputlevel]} {set TAMDoutputlevel 3}
