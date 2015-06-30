# tcl library for CFACV
# Cameron F Abrams
# 2009-14

# NAMD requeriments:
#
# puts must be replaced by print. Other commands
# could be found in NAMD user guide.

proc cfacv_banner { argv } {
    global CFACV_VERSION
    global CFACV_BASEDIR
    print "================================================"
    print "CFACV) Collective Variables in VMD/NAMD v. $CFACV_VERSION"
    cfacvBanner
    print "CFACV) 2009-2011, Cameron F Abrams, Drexel University"
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
###########################################################
# Primitive data handling procedures
###########################################################

# ListToArray: allocates a new array and assigns it values
# from the Tcl list $l; returns pointer to new array.  List
# elements are treated as double-precision by default.
#
proc ListToArray {l} {
    set length [llength $l]
    set a [new_array $length]
    set i 0
    foreach item $l {
        array_setitem $a $i $item
        incr i 1
    }
    return $a
}

# ListToArray_Data:  Assigns elements of an existing 
# array pointed to by $a from the tcl list $l
#
proc ListToArray_Data { a l } {
    set i 0
    foreach item $l {
	array_setitem $a $i $item
	incr i 1
    }
    return $a
}

# intListToArray: like listToArray except elements of
# $l are treated as integers
#
proc intListToArray {l} {
   set length [llength $l]
   set a [new_arrayint $length]
   set i 0
   foreach item $l {
        arrayint_setitem $a $i $item
        incr i 1
    }
  return $a
}

# intListToArray_Data: list listToArray_Data except elements
# of $l are treated as integers
#
proc intListToArray_Data {a l} {
   set length [llength $l]
   set i 0
   foreach item $l {
        arrayint_setitem $a $i $item
        incr i 1
    }
  return $a
}

# ArrayToList:  Creates a new tcl list of length $n and
# assigns it elements of existing array $a
#
proc ArrayToList {a n} {
    set length $n
    set l {}
    for {set i 0} { $i < $length} {incr i} {
	lappend l [array_getitem $a $i]
    }
    return $l
}

# intArrayToList: like ArrayToList but treats
# elements of $a as integers
#
proc intArrayToList {a n} {
    set length $n
    set l {}
    for {set i 0} { $i < $length} {incr i} {
	lappend l [arrayint_getitem $a $i]
    }
    return $l
}


###############################################################
# read_centersPDB
###############################################################
# Sets up the centers according to the template PDB file
# in which each atom's beta field designates the index of
# the center to which it belongs; each atom can belong
# to one and only one center!  Also, the occupancy
# field should contain the atomic mass in amu.
#
# IMPORTANT: A beta of i means the atom belongs to the (i-1)th 
# center! A beta of "0" indicates that the atom does not belong 
# to a center.
# 
# returns:
#  serArr: a tcl list in which the i'th element
#          is a list of atom serial numbers (1-based indicies)
#          belonging to center i
#
#          e.g.
#          {1 2 4    } ;#group 1
#          {3 5 7    } ;#group 2
#          {6 9 12 11} ;#group 3
#
#  masses: a tcl list in which the i'th element
#          is the total mass of center i
#
#          e.g.
#          {1.004 15.99 1.004       } #group 1
#          {1.004 15.99 1.004       } #group 2
#          {1.004 15.99 1.004 1.004 } #group 3
#
#  return value is the number of centers
#
###############################################################
proc read_centersPDB { templatePdb serArr mass pdb} {
    upvar $serArr serArray
    upvar $mass masses
    upvar $pdb pdbline 

    # pdbline is a list with all the properties from the pdb file
    # this allow to access to the resname or segname later, which will
    # allow to use resname or segname keywords in the cv declaration

    # Open and read the pdb. Put all it in lines variable.
    set inStream [open $templatePdb r]
    set data [read $inStream]
    close $inStream
    set lines [split $data \n]

    # Get list indices with the center numbers and Nmon as the number of
    # centers
    set indices {}
    foreach line $lines {
	set cardstr [string range $line 0 3]
	set betastr [string trim [string range $line 60 65]]
	# if the beta field is nonzero, it is assumed to be an integer
	# that indicates which center this atom belongs to
	if {([string equal $cardstr "ATOM"]) && [expr 0==[string equal $betastr "0.00"]]} {
          lappend indices [expr int($betastr)-1]
        }
    }
    set indices [lsort -integer -unique $indices]
    set nMon [llength $indices]

    # contruct the void list of masses, atom serial and reslist
    for {set i 0} {$i < $nMon} {incr i} {
	lappend serArray {}
	lappend pdbline {}
	lappend masses 0.0
        lappend reslist 0
    }

    # append list of masses and serials ID to each group index
    set n 0
    foreach line $lines {
	set cardstr [string range $line 0 3]
	set serlstr [string trim [string range $line 5 11]]
	set resname [string trim [string range $line 17 20]]
	set occustr [string trim [string range $line 54 59]] ;# used to hold atomic mass in amu.
	set betastr [string trim [string range $line 60 65]]
	if {([string equal $cardstr "ATOM"]) && [expr 0==[string equal $betastr "0.00"]]} {
	    
	    set t [expr int($betastr)-1] ; #Center index
	    set nlist [lindex $serArray $t] ; # get atom index list for center $t

	    # add this index to the list; index is serial-1
	    # since begining with 100000, serial numbers are stored and
	    # read in hex, need to check whether this is a hex number
	    # -- I'll do this by seeing if it is not an integer
#	    print "adding serl $serlstr to $nlist"
	    if {[string is integer $serlstr]} {
		set serl [expr $serlstr]
	    } else {
		set serl [expr 0x$serlstr]
	    }
	    lappend nlist $serl
	    # return the newly appended list of indices to the list of lists
	    set serArray [lreplace $serArray $t $t $nlist]

            # update all the properties. Now, if is a group, pdbline
            # correspond to the last line of the group. TODO I don't know it
            # have sense for a group this mechanism, but it is good to
            # select many atoms in the cv.  I guess that it should be better
            # to put one entrance per atom in the pdbline, then the info of
            # the groups will be saved, and avoid serarray and any other
            # stuff that complicate the legibility
	    set pdbline [lreplace $pdbline $t $t "$cardstr $serlstr $resname $occustr $betastr"]

	    # update the mass of this pseudoatom
	    set am [expr double($occustr)]
	    set masses [lreplace $masses $t $t [expr [lindex $masses $t] + $am]]

            # Renames
            set reslist [lreplace $reslist $t $t $resname]

	}
	incr n
    }

    # FIXME: This code was here to allow compute chapeau functions separatedly
    # for different pair types of particles. For instance, this allow to
    # recover SOD SOD, CLA CLA and SOD CLA pair potentials in 1 TAMD
    # simulation. Each index has a number in ch_id which allow to sort the pair
    # in the different chapeau objects on the c code.  From the studies with
    # SOD CLA, this pair potentials will be OK only if the ficticius
    # temperature is the same that the real one.  On the other hand, a better
    # way to achive this is needed (without saving a lot of numbers in ch_id).
    # For understand how this worked, see the previous versions of the code.
    # set k 0
    # set test [expr 0==[string equal $chlist "{}"]]
    # for {set i 0} {$i < $nMon} {incr i} {
    #   for {set j [expr $i+1]} {$j < $nMon} {incr j} {

    #     if $test {
    #       set pair1 "[lindex $reslist $i] [lindex $reslist $j]" ;# e.g. "SOD CLA"
    #       set pair2 "[lindex $reslist $j] [lindex $reslist $i]" ;# e.g. "CLA SOD"
    #       set a [lsearch $chlist $pair1]
    #       set b [lsearch $chlist $pair2]
    #       
    #       set k [expr {$a>$b? $a: $b}]
    #     }

    #     # with the next expresion the [expr $j+($nCntr-2)*$i-($i-1)*$i/2-1] element of ch_id is
    #     # -1 for non intereset pair, or the corresponding n-index or the $chlist.
    #     lappend ch_id $k
    #   }
    # }

    # if {([string equal [lsort -unique $ch_id] "-1"])} {
    #   error "ERROR: No residues found to accelerate"
    # }

    print "DB: returning nMon $nMon"
    return $nMon
}


proc center_selections { molID serArray } {
    set p {}
    set n [llength $serArray]
    print "DB: center_selections thinks there are $n centers"
    set i 0
    for {set i 0} { $i < $n } { incr i } {
	set arr [lindex $serArray $i]
	set a [atomselect $molID "serial $arr"]
        $a global
	lappend p $a
    }
    return $p
}

######################################################################
# read_centersVMD
#
# If a molecule is loaded in VMD and a "centers.dat" file
# (generated by mk_tPDB.tcl) exists, then it is possible
# to define the atomselections for each center without
# reading a template PDB.  Because "atomselections" are
# not part of the Tcl interface to NAMD, this approach
# is not useful for tcl_forces.  
#
# pL is a Tcl list of atomselections; this procedure only appends it.
# fileName is the name of the "centers.dat" file
# molID is the molecule id number of the molecule for
# which centers are defined.
#
# The format of the "centers.dat" file is such that
# one center is defined on each line.  The format
# of a line is
#
# id mass residue_1 residue_2 ... residue_N
#
# where id is the index of the center, mass is its mass in amu
# and residue_i is a string concatenation of the single-letter
# amino acid designation and the residue number.
#
#####################################################################
proc read_centersVMD { pL mL fileName molID } {
    upvar $pL p
    upvar $mL m
    set ofp [open "$fileName" "r"]
    set data [read -nonewline $ofp]
    close $ofp
    set lines [split $data \n]
    foreach l $lines {
	set il {}
	lappend m [lindex $l 1]
	set RL [lreplace $l 0 1]
	foreach ri $RL {
	    lappend il [string range $ri 1 end]
	}
#	print "DB: making selection from residues $il"
	set a [atomselect $molID "protein and resid $il"]
	$a global
	lappend p $a 
    }
    return [llength $lines]
}

proc read_centers_residLists { ridL fileName } {
    upvar $ridL r
    set ofp [open "$fileName" "r"]
    set data [read -nonewline $ofp]
    close $ofp
    set lines [split $data \n]
    foreach l $lines {
	set il {}
	set RL [lreplace $l 0 1]
	foreach ri $RL {
	    lappend il [string range $ri 1 end]
	}
	lappend r $il 
    }
    return [llength $lines]    
}

proc read_cvs { cv_file cv_list pdb} {
  #and then cv_list is the append of dat.... why all the stuff?
  upvar $pdb pdbline 
  upvar $cv_list cvL

  set cvL {}

  # read the cv file
  set inStream [open $cv_file "r"]
  set data [read -nonewline $inStream]
  close $inStream
  set lines [split $data \n]

  set ncv 0
  foreach line $lines {

    # Allow empty lines and commentries
    if {[string index $line 0] == "#"} continue
    if {[string length $line] < 0} continue

    # Supouse that the line is "CARTESIAN_X   0"...
    set dat {}
    foreach l $line { lappend dat $l } ;#... dat is "CARTESIAN_X 0"
    set typ [lindex $dat 0]            ;#... typ is "CARTESIAN_X"
    set ind [lreplace $dat 0 0]        ;#... ind is "0"

    # Search for optional strings after the index list
    set optlist {}
    for {set i 0} {[string is integer [lindex $ind $i]] && $i < [llength $ind]} { incr i } {}
    if {$i < [llength $ind]} {
        # Put in optlist that string and let in ind only the numbers
        set optlist [lreplace $ind 0 [expr $i-1]]
        set ind [lreplace $ind $i end]
    }
    set nInd [llength $ind]
    set nOpt [llength $optlist]
    
    # Si no hay indices, busco palabras clave de seleccion
    if {$nInd==0} {
      set j 0
      foreach o $optlist {
        incr j

        # fill index by resname
        if {$o=="resname"} {
          print [lindex $optlist  $j]
          set resn [lindex $optlist  $j]
          foreach l $pdbline {
            set index [lindex $l 4]
            if {$resn==[lindex $l 2]} {lappend ind [expr $index-1]}
          }

          # Remove resname from option list
          set j [expr $j-1]
          set optlist [lreplace $optlist $j $j ]
          set optlist [lreplace $optlist $j $j ]

          break
        }

      }
    }
    set nInd [llength $ind]

    # Si todavia no hay indices, doy error
    if {$nInd==0} {error "Error, no indices para la linea $line"}


    lappend cvL [list $typ $ind $optlist] ;#... cvL is cvL+dat
    incr ncv
    
    # Check for a valid CV
    switch $typ {
      ZSDCIRCLE {
        set_zsd_circle [lindex $optlist 0] [lindex $optlist 1]  [lindex $optlist 2] [lindex $optlist 3]
      }
      ZSDXRANGE {
        set_zsd_circle [lindex $optlist 0] 0.0  [lindex $optlist 1] [lindex $optlist 2]
      }
      ZSDRING {
        set_zsd_ring [lindex $optlist 0] [lindex $optlist 1]  [lindex $optlist 2] [lindex $optlist 3] [lindex $optlist 4] 
      }
      BOND -
      S -
      ANGLE -
      DIHED -
      COGX -
      COGY -
      COGZ -
      CARTESIAN_X -
      CARTESIAN_Y -
      CARTESIAN_Z {}
      default {error "ERROR: $typ is not a valid CV type in $cv_file."}
    }
  }

  return $ncv
}

proc generate_cv { cvOpt p } {

    puts -nonewline "CFACV) Generating cv.inp from $cvOpt..."; flush stdout
    set retList {}
    set segidList {}

    set nCntr [llength $p]
    
    foreach s $p {
	set segid [lsort -unique [$s get segid]]
	if {[llength $segid] > 1} {
	    print "ERROR: a center owns atoms from more than one segid: $segid"
	}
	lappend segidList $segid
    }

#    print "DB: generate_cv: segidList: $segidList" 

    if {($cvOpt == "autogenerate-bonds" || $cvOpt == "all-to-all-bonds")} {
	for {set ii 0} {$ii < [expr $nCntr-1]} { incr ii } {
	    for {set jj [expr $ii + 1]} { $jj < $nCntr } { incr jj } {
		set thisind [list $ii $jj]
		lappend retList [list "BOND" $thisind]
	    }
	}
    } elseif {[string first "shortest" $cvOpt]!=-1} {
	set tl [split $cvOpt -]
	set nIntraSeg [lindex $tl 1]
	if {[llength $tl] > 2} {
	    set nInterSeg [lindex $tl 2]
	}
	set ii 0
	set jj 0
	set intraSegL {}
	set interSegL {}
	foreach s1 $p {
	    set iseg [lindex $segidList $ii]
	    set c1 [measure center $s1 weight mass]
	    set jj [expr $ii + 1]
	    foreach s2 [lreplace $p 0 $ii] {
		set jseg [lindex $segidList $jj]
		set c2 [measure center $s2 weight mass]
		set d12 [veclength [vecsub $c1 $c2]]
		if {$iseg == $jseg} {
		    lappend intraSegL [list $ii $jj $d12]
		} else {
		    lappend interSegL [list $ii $jj $d12]
		}
		incr jj
	    }
	    incr ii
	}
	set intraSegL [lsort -real -index 2 $intraSegL]
	set interSegL [lsort -real -index 2 $interSegL]
	for {set i 0} {$i < $nIntraSeg} {incr i} {
	    set tl [lindex $intraSegL $i]
	    set ind [list [lindex $tl 0] [lindex $tl 1]]
	    lappend retList [list "BOND" $ind]
	}
	if {[info exists nInterSeg]} {
	    for {set i 0} {$i < $nInterSeg} {incr i} {
		set tl [lindex $interSegL $i]
		set ind [list [lindex $tl 0] [lindex $tl 1]]
		lappend retList [list "BOND" $ind]
	    }
	}
	set retList [lsort -index 1 $retList]
    } elseif {[string first "maxlength" $cvOpt] != -1} {
	set tl [split $cvOpt =]
	set maxL_intra [lindex $tl 1]
	if {[llength $tl] > 2} {
	    set maxL_inter [lindex $tl 2]
	}
	set ii 0
	set jj 0
	set intraSegL {}
	set interSegL {}
	foreach s1 $p {
	    set iseg [lindex $segidList $ii]
	    set c1 [measure center $s1 weight mass]
	    set jj [expr $ii + 1]
	    foreach s2 [lreplace $p 0 $ii] {
		set jseg [lindex $segidList $jj]
		set c2 [measure center $s2 weight mass]
		set d12 [veclength [vecsub $c1 $c2]]
		if {$iseg == $jseg} {
		    if {[expr $d12 < $maxL_intra]} {
			lappend intraSegL [list $ii $jj $d12]
		    }
		} else {
		    if {[expr $d12 < $maxL_inter]} {
			lappend interSegL [list $ii $jj $d12]
		    }
		}
		incr jj
	    }
	    incr ii
	}
	set intraSegL [lsort -real -index 2 $intraSegL]
	set interSegL [lsort -real -index 2 $interSegL]
	set nIntraSeg [llength $intraSegL]
	set nInterSeg [llength $interSegL]
	for {set i 0} {$i < $nIntraSeg} {incr i} {
	    set tl [lindex $intraSegL $i]
	    set ind [list [lindex $tl 0] [lindex $tl 1]]
	    lappend retList [list "BOND" $ind]
	}
	if {[info exists nInterSeg]} {
	    for {set i 0} {$i < $nInterSeg} {incr i} {
		set tl [lindex $interSegL $i]
		set ind [list [lindex $tl 0] [lindex $tl 1]]
		lappend retList [list "BOND" $ind]
	    }
	}
	set retList [lsort -index 1 $retList]
    } elseif {[string first "cartesian" $cvOpt] != -1} {
	for {set ii 0} {$ii < $nCntr} { incr ii } {
	    foreach c {X Y Z} {
		lappend retList [list "CARTESIAN_${c}" $ii]
	    }
	}
    } elseif {[string first "peptide_tilt" $cvOpt] != -1} {
	set nDihed [expr ($nCntr/2)-1]
	# O 1, 3, 5, 7, 9, etc.
        # H 2, 4, 6, 8, 10, etc.
	# dihed: 2-1-3-4, 4-3-5-6, 6-5-7-8, 8-7-9-12, etc.
	# but begin counting at 0!!!
	set p1 1
	set p2 0
	set p3 2
	set p4 3
	for {set i 0} {$i < $nDihed} {incr i} {
	    lappend retList [list "DIHED" $p1 $p2 $p3 $p4]
	    incr p1 2
	    incr p2 2
	    incr p3 2
	    incr p4 2
	}
    } elseif {[string first "phi-psi" $cvOpt] != -1} {
	set nPsi [expr ($nCntr/3)-1]
	print "CFACV) [expr ($nCntr/3)] residues -> $nPsi psi angles"
	set psia 0
	set psib 1
	set psic 2
	set psid 3
	set phia 2
	set phib 3
	set phic 4
	set phid 5
	for {set i 0} {$i < $nPsi} {incr i} {
	    lappend retList [list "DIHED" $psia $psib $psic $psid]
	    lappend retList [list "DIHED" $phia $phib $phic $phid]
	    incr psia 3
	    incr psib 3
	    incr psic 3
	    incr psid 3
	    incr phia 3
	    incr phib 3
	    incr phic 3
	    incr phid 3
	}
    } else {
	print "ERROR: cvOpt $cvOpt is not recognized."
    }
    print "done."
    return $retList
}

proc output_cv { cvList fileName } {
    set fp [open "$fileName" "w"]
    foreach cv $cvList {
	set typ [lindex $cv 0]
	set ind [lreplace $cv 0 0]
	print $fp "$typ $ind"
    }
    close $fp
}

####################################################################
# Procedures for handling restraints
####################################################################
#
# A "restraint" is the value of a linear combination of collective
# variables:
# 
# R_l = \sum_{k=1}^{nCV} C_{kl} \Theta_k
#
# where \Theta_k is the value of the k'th collective variable, and 
# C_{kl} is the coefficient multipliying \Theta_k in the l'th
# restraint.
#
# A restraint measure is 
# 
# Z_l = R_l - B_1
# 
# where B_1 is the target value of the restraint.
#
# The potential energy of the restraint is harmonic:
#
# U_l = (1/2) k_l Z_l^2
#
# so the force on a degree of freedom x is
#
# f = -\grad_x U_l = -k_l Z_l \grad_x Z_l
#
#
# Example:  If collective variable \Theta_0 is a bond between centers 0 and 1
# that we wish to restrain to a distance of B_0 = 10 A using restraint R_0, and there
# are a total of 3 collective variables defined in the system, then
# R_0 = (1)(\Theta_0) + (0)(\Theta_1} + (0)(\Theta_2)
# and
# Z_0 = R_0 - B_0
#
# This allows one to define restraints that are linear combinations
# of collective variables.  For cases in which only isolated collective variables
# are restrained, the matrix C_{k1} is the NxN identity matrix where N is
# the number of collective variables.
# 
proc cvc_getcvi { cvc } {
    set r {}
    set i 0
    foreach cv $cvc {
	if {[expr $cv > 0.0]} {
	    lappend r $i
	}
	incr i
    }
    return $r
}



proc rlist_setup { ncv } {
  global ds
  global restr
  global temperature
  global TAMDverbose

  for {set i 1} {$i<=$restr(num)} {incr i} {

    # Set cvc
    set restr($i.cvc) {}
    for {set j 0} { $j < $ncv} {incr j} {lappend restr($i.cvc) 0}
    foreach cvindex [lsort -unique $restr($i.cv)] {
      set j [expr $cvindex-1]
      set restr($i.cvc) [lreplace $restr($i.cvc) $j $j 1]
    }


    # Defaults and consistency in function
    if {![info exists restr($i.func)]} {set restr($i.func) HARMONIC}
    set restr($i.func) [string toupper $restr($i.func)]
    if {![info exists restr($i.k)]} {print "CFACV) WARNING: restraint with spring constant 0"; set restr($i.k) 0}


    # Defaults and consistency in boundaries
    if {![info exists restr($i.bound)]} {set restr($i.bound) NADA}
    set restr($i.bound) [string toupper $restr($i.bound)]
    switch $restr($i.bound) {
      NADA {}
      PERIODIC {
        set restr($i.min) [expr -1*acos(-1)]
        set restr($i.max) [expr acos(-1)]
      }
      PBC {
        if {![info exists restr($i.min)]}    {error "CFACV) ERROR: restraint need max value"}
        if {![info exists restr($i.max)]}    {error "CFACV) ERROR: restraint need min value"}
      }
      SOFTUPPER - SOFTLOWER - SOFT {
        if {![info exists restr($i.min)]}    {error "CFACV) ERROR: restraint need max value"}
        if {![info exists restr($i.max)]}    {error "CFACV) ERROR: restraint need min value"}
        if {![info exists restr($i.boundk)]} {error "CFACV) ERROR: restraint need boundk value"}
      }
      default {error "ERROR: $restr($i.bound) is not a valid boundary type."} 
    }


    # Allocating ds->restr[i] in the C code.
    set pcvc [ListToArray $restr($i.cvc)]
    set inicial 0 ; # obsolet when r->z=r->val is placed at first iteration (see cfacv.c)
    set restr($i.address) [DataSpace_AddRestr $ds $restr($i.k) $inicial $ncv $pcvc $restr($i.func) $restr($i.min) $restr($i.max) $restr($i.bound) $restr($i.boundk)]
     
    # Default and consistency in type. Allocating ds->restr[i]->tamdOpt in the C code.
    if {![info exists restr($i.type)]} {error "CFACV) ERROR: restraint need a type"}
    set restr($i.type) [string toupper $restr($i.type)]
    switch $restr($i.type) {
      TAMD {
        if {![info exists restr($i.g)]}     {error "CFACV) ERROR: TAMD restraint need g"}
        if {![info exists restr($i.dt)]}    {error "CFACV) ERROR: TAMD restraint need dt"}
        if {![info exists restr($i.temp)]}  {error "CFACV) ERROR: TAMD temperature variable is needed."}
        if {![info exists restr($i.chid)]}  {set restr($i.chid) 0}

        # DataSpace_AddTamdOpt $ds $ir $restr($i.g) [expr $temperature*$restr($i.biasf)*0.00198735] $restr($i.dt)
        restr_AddTamdOpt $restr($i.address) $restr($i.g) [expr $restr($i.temp)*0.00198735] $restr($i.dt)
      }
      SMD {
        if {![info exists restr($i.target)]} {error "CFACV) ERROR: SMD restraint need target"}
        if {![info exists restr($i.ini)]} {error "CFACV) ERROR: SMD restraint need ini"}
        if {![info exists restr($i.end)]} {error "CFACV) ERROR: SMD restraint need end"}
        restr_AddSmdOpt $restr($i.address) $restr($i.target) $restr($i.ini) $restr($i.end)
      }
      default {error "ERROR: $restr($i.type) is not a valid boundary type."} 
    } 

  }

  if {[info exists TAMDverbose]} {foreach ri [array get restr] {print "CFACV) $ri"}}
  
}



proc my_getcellsize { XSCFILE } {
    if {[info exists XSCFILE] && [file exists $XSCFILE]} {
    } else {
      error "CFACV) ERROR: no xsc file for box info."
    }
    #print "CFACV) DEBUG my_getcellsize opening $XSCFILE"
    set fp [open $XSCFILE "r"]
    set dat [read -nonewline $fp]
    #print "CFACV) DEBUG my_getcellsize dat:"
    #print "$dat"
    set lines [split $dat \n]
    #print "CFACV) DEBUG my_getcellsize lines:"
    #print "$lines"
    set lastline [lindex $lines end]
    lappend L [lindex $lastline 1]
    lappend L [lindex $lastline 5]
    lappend L [lindex $lastline 9]
    print "CFACV) INFO Box size $L"
    close $fp
    return $L
}

proc my_getorigin { XSCFILE } {
    if {[info exists XSCFILE] && [file exists $XSCFILE]} {
    } else {
      error "CFACV) ERROR: no xsc file for box info."
    }
    set fp [open $XSCFILE "r"]
    set dat [read -nonewline $fp]
    set lines [split $dat \n]
    set lastline [lindex $lines end]
    lappend O [lindex $lastline 10]
    lappend O [lindex $lastline 11]
    lappend O [lindex $lastline 12]
    print "CFACV) INFO origin $O"
    close $fp
    return $O
}

proc assign_center_betas { p includeH } {
    set nCntr [llength $p]
    for {set i 0} { $i < $nCntr } { incr i } {
	set cntrIndex [expr $i + 1]
	set a [lindex $p $i]
	if { $includeH } {
	    $a set beta $cntrIndex
	} else {
	    set thisbeta {}
	    foreach e [$a get element] {
		if { $e != "H" } {
		    lappend thisbeta $cntrIndex
		} else {
		    lappend thisbeta 0.0
		}
	    }
	    $a set beta $thisbeta
	}
	$a set occupancy [$a get mass]
    }
}

proc report_centers { p fileName } {
    set nCntr [llength $p]
    set ofp [open "$fileName" "w"]
    for {set i 0} { $i < $nCntr } { incr i } {
	set a [lindex $p $i]
	set this_m [format "%8.2f" [vecsum [$a get mass]]]
	puts -nonewline $ofp "[format %3i $i] $this_m "
	foreach id [lsort -integer -unique -index 1 [$a get {resname resid chain}]] {
	    set indx "[aa_321 [lindex $id 0]][format %-3i [lindex $id 1]]"
	    puts -nonewline $ofp "$indx "
	}
	puts $ofp ""
    }
    close $ofp
}

proc initialize_base_template { pdb } {
    mol new  $pdb

    set molid [molinfo top get id]

    set all [atomselect top all]
    $all set beta 0
    $all set occupancy 1.0

    $all delete
   
    return $molid
}

proc within_range {selList selB r} {
  set nc 0
  foreach sel $selList {
      set ABlist [measure contacts $r $sel $selB]
      set Alist [lindex $ABlist 0]
      set Blist [lindex $ABlist 1]
      if {[llength $Alist] > 0} {
          incr nc [llength $Alist]
      }
  }
  return $nc
}
 
proc cvs_setup { } {
  global cvList
  global ds

  foreach cv $cvList {
    set typ [lindex $cv 0]
    set pind [intListToArray [lindex $cv 1]]
    set nind [llength [lindex $cv 1]]

    # Assuming cv=CARTESIAN_X 0 then typ=CARTESIAN_X, pind={0}, nind=1
    DataSpace_AddCV $ds $typ $nind $pind 
  }

}

proc chapeau_setup { NUMREP ds XSCFILE min N max begin_evolve usetamdforces outputfile outfreq outlevel nupdate} {
    global chape

    # FIXME: This code was here to allow compute chapeau functions separatedly
    # for different pair types of particles. For instance, this allow to
    # recover SOD SOD, CLA CLA and SOD CLA pair potentials in 1 TAMD
    # simulation. Each index has a number in ch_id which allow to sort the pair
    # in the different chapeau objects on the c code.  From the studies with
    # SOD CLA, this pair potentials will be OK only if the ficticius
    # temperature is the same that the real one.  On the other hand, a better
    # way to achive this is needed (without saving a lot of numbers in ch_id).
    # For understand how this worked, see the previous versions of the code.
    #
    # # chplist holds the list of pair potential types 
    # # e.g. {{SOD CLA} {CLA CLA} {SOD SOD}}
    # if {![info exists chlist]} {set chlist "{}"}
    # set chnum [llength $chlist]
    #
    # # Pair potentinal interaction
    # intListToArray_Data [DataSpace_chid $ds ] $ch_id
    # print "CFACV) The pair chapeau (map in vector) is:"
    # for {set i 0} { $i < $nCntr } { incr i } {
    #   set aux ""
    #   for {set j 0} { $j <= $i } { incr j } {set aux "$aux  "}
    #   for {set j [expr $i+1]} { $j < $nCntr } { incr j } {
    #     set aux "$aux [lindex $ch_id [expr $j+($nCntr-2)*$i-($i-1)*$i/2-1]]"
    #   }
    #   puts $aux
    # }

    # FIXME. Add a index to load a  initial knots file for each chapeau
    # if {[info exists initKnotsINP]} {Tcl_DataSpace_InitKnots $ds $initKnotsINP}

    # FIXME: pairwise analytical potential
    # DataSpace_SetupPairCalc $ds $min $N $coff $nlcoff $begin_evolve $usetamdforces $outputfile $outfreq $outlevel $nupdate $chnum

    # TODO: I guess one setup subroutine per chapeau object would be nice, so each chapeau can have different grid
    # and be asociated with differen CV sets

    # Currently only option is a 1D grid
    DataSpace_Setup1Dchapeau $ds $NUMREP $min $N $max $begin_evolve $usetamdforces $outputfile $outfreq $outlevel $nupdate

    # Get the adress of the chapeaus
    for {set i 0} {$i<$NUMREP} {incr i} {
      set chape($i.address) [DataSpace_get_chapeauadress $ds $i]
    }

    if [string equal $XSCFILE "off"] {
      print "CFACV) DEBUG: Tcl_InitializePairCalc nobox"
      DataSpace_SetupPBC $ds 0  0 0 0  0 0 0
    } else {
      set LL [my_getcellsize $XSCFILE]
      set O [my_getorigin $XSCFILE]
      print "CFACV) DEBUG: Tcl_InitializePairCalc box size [lindex $LL 0] [lindex $LL 1] [lindex $LL 2]"
      print "CFACV) DEBUG: Tcl_InitializePairCalc origin   [lindex $O 0] [lindex $O 1] [lindex $O 2]"
      DataSpace_SetupPBC $ds 1 [lindex $O 0] [lindex $O 1] [lindex $O 2] [lindex $LL 0] [lindex $LL 1] [lindex $LL 2] 
    }
}

proc Tcl_Reinitialize { ds restartINP } {
    if {[file exists $restartINP]} {
	print "CFACV) Reinitializing Z values from $restartINP"
	set fp [open $restartINP "r"]
	set zinp [read -nonewline $fp]
	close $fp
	DataSpace_SetRestraints $ds [ListToArray $zinp]
	return 0
    } else {
	print "CFACV) Error: could not find restart file $restartINP"
	print "CFACV) Initial z-values taken from CV's."
	return 1
    }
}

proc Tcl_UpdateDataSpace { ds lC groups first timestep } {
    upvar $lC p

    # Move group center positions to dataspace
    set i 0
    foreach g $groups {
	ListToArray_Data [DataSpace_centerPos $ds $i] $p($g)
	incr i
    }
    MyParanoiaCheck $ds "tripped after moving data to dataspace"

    # Compute CV's within dataspace
    DataSpace_ComputeCVs $ds
    MyParanoiaCheck $ds "tripped after computing CV's"
    # At this point, the (x,y,z) position data is no longer needed.
    # We can now write into its space the (x,y,z) restraint forces.

    # Compute restraining forces
    DataSpace_RestrainingForces $ds $first $timestep    
    MyParanoiaCheck $ds "tripped after computing restraining forces"

    #  fprintf(stderr,"AAAAAAAAAAAS");

    # Transmit forces back to groups
    set i 0
    foreach g $groups {
        addforce $g [ArrayToList [DataSpace_centerPos $ds $i] 3]
        incr i
    }

    # Add restraint energy to NAMD energy structure
    addenergy [DataSpace_RestraintEnergy $ds]
}

proc Tcl_ObserveDataSpace { ds cntrSel frame } {
    set i 0
    foreach sel $cntrSel {
	$sel frame $frame
#	print "DB: Tcl_ObserveDataSpace center $i [measure center $sel weight mass]"
	ListToArray_Data [DataSpace_centerPos $ds $i] [measure center $sel weight mass]
	incr i
    }
    DataSpace_ComputeCVs $ds
}

# This routine uses the DataSpace_checkdata function
# to detect any "infs" or "nans" in the position data
# and dies with output if any are detected.
# Note this is a stub if a variable name PARANOIA
# is *not* defined in either the NAMD config file
# or the cfacv_tclforces script.
proc MyParanoiaCheck {ds msg} {
    global PARANOIA
    if {[info exists PARANOIA]} {
	if {$PARANOIA} {
	    if {[DataSpace_checkdata $ds]} {
		print "CFACV/PARANOIA) $msg"
		DataSpace_dump $ds
		exit
	    }
	}
    }
}

proc Tcl_DataSpace_InitKnots { ds filename } {
    if {[file exists $filename]} {
	DataSpace_InitKnots $ds $filename
    } {
	error "CFACV) error: initial knots file $filename not found."
    }

}
