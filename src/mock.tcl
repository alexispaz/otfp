# This Tcl script provides a **mock of NAMD** for testing and development purposes.
#
# The script simulates a typical NAMD workflow loading a coordinate file,
# runing OTFP to compute extra forces, and advancing frames, without requiring
# the actual NAMD executable.
#
# Usage:
#   tclsh mock.tcl <xyz_base> <namd_script>
#
#   where:
#     <xyz_base>    = base name of the coordinate file (e.g., job0)
#     <namd_script> = NAMD script to simulate (e.g., OTFP.namd)
#
# Example:
#   tclsh ../../src/mock.tcl job0 OTFP.namd
#
# Copyright (c) 2025 Sergio Alexis Paz <apaz@unc.edu.ar>

namespace eval ::mock {
    
  # Result from OTFP code
  # results(g0.force) is the force on group g0 
  # results(energy) is the energy to be added
  variable results
  array unset results

  # Info from previous NAMD run. For instance:
  # atoms(g0) are the atom indices of group g0 
  # atoms(a1.sym) are the chemical symbol of atom 1
  # atoms(na) are the number of atoms
  # atoms(f1.a1.pos) are the position of atom 1 in frame 1
  # atoms(nf) are the numbers of frames
  variable atoms 
  array unset atoms

  # Group partition from addgroup
  variable groups
  set groups {} 

  # tclforcesscript
  variable script

  # outputname and mock files
  variable oname
  variable f1
  variable f2

  # structure  
  variable pdbfile

  variable step
  set step 0
         
  # Stubs for process used in OTFP source
  proc tclforces args {}
  proc enabletotalforces args {}
  proc loadforces args {}
  proc loadtotalforces args {}

  # Automatic stubs of unused namd commands from input name
  proc make_stubs {filename} {

      # Detectar el namespace desde donde se llamó make_stubs
      set caller_ns [uplevel 1 {namespace current}]

      # Leemos el contenido
      set f [open $filename r]
      set content [read $f]
      close $f

      # Obtenemos todas las "palabras" que aparecen al inicio de línea (posibles comandos)
      foreach line [split $content "\n"] {
          if {[regexp {^\s*#} $line]} { continue } ;# saltar comentarios
          if {[regexp {^\s*$} $line]} { continue } ;# saltar líneas vacías

          # Extraer el primer token (comando)
          if {[regexp {^\s*([A-Za-z0-9_-]+)} $line -> cmd]} {

              # Si el comando no existe en Tcl, lo definimos como stub vacío
              set fullcmd ${caller_ns}::$cmd
              if {![llength [info commands $cmd]]} {
                  puts "Definiendo stub para $fullcmd"
                  proc $fullcmd {args} {}
              }
          }
      }

  }

  # Mock of NAMD commands
     
  proc coordinates {filename} {
    variable pdbfile
    variable atoms

    set pdbfile $filename
    set f [open $filename r]
    set i 1
    while {[gets $f line] >= 0} {
        if {[string match "ATOM*" $line] || [string match "HETATM*" $line]} {
            # Leer los campos del formato PDB estándar
            set atom_name [string trim [string range $line 12 15]]
            set res_name  [string trim [string range $line 17 19]]
            set chain_id  [string trim [string range $line 21 21]]
            set res_id    [string trim [string range $line 22 25]]
            set x [string trim [string range $line 30 37]]
            set y [string trim [string range $line 38 45]]
            set z [string trim [string range $line 46 53]]

            # Guardar los datos
            set atoms(a$i.sym) $atom_name
            set atoms(a$i.resname) $res_name
            set atoms(a$i.chain) $chain_id
            set atoms(a$i.resid) $res_id

            set atoms(f0.a$i.pos) [list $x $y $z]

            incr i
        }
    }
    close $f
    set atoms(na) $i

  }

  proc binCoordinates {filename} {
    variable pdbfile
    exec catdcd -o "mock.tmp" -s $pdbfile -otype pdb -namdbin $filename
    coordinates "mock.tmp"
  }

  proc tclforcesscript {scriptfile} {
      variable script
      set script $scriptfile
  }

  proc print {thing} {
      puts $thing
  }

  proc getstep {} {
      variable step
      return $step
  }

  proc addgroup {atom_indices} {
      variable groups
      variable atoms

      set gid "g[llength $groups]"
      set atoms($gid) $atom_indices
      lappend groups $gid
      return $gid
  }


  proc loadcoords {arrName} {
      variable groups
      variable atoms
      variable step

      upvar $arrName arr
      array unset arr 2>/dev/null

      foreach g $groups {
          set sumx 0.0
          set sumy 0.0
          set sumz 0.0

          foreach a $atoms($g) {
              set x [lindex $atoms(f$step.a$a.pos) 0]
              set y [lindex $atoms(f$step.a$a.pos) 1]
              set z [lindex $atoms(f$step.a$a.pos) 2]
              set sumx [expr {$sumx + $x}]
              set sumy [expr {$sumy + $y}]
              set sumz [expr {$sumz + $z}]
          }

          # centro de masa
          set n [llength $atoms($g)]
          set com [list [expr {$sumx/$n}] [expr {$sumy/$n}] [expr {$sumz/$n}]]
          set arr($g) $com
      }
  }


  proc addforce {gid vec} {
      variable results
      set results($gid.force) $vec
  }
     
  proc addenergy {e} {
      variable results
      set results(energy) $e
  }  


  proc outputName {filename} {
    variable oname
    variable f1
    variable f2
    set oneame $filename

    loadxyz "$::oname.xyz"
             
    set f1 [open "$::oname.mk1" w]
    set f2 [open "$::oname.mk2" w]

  }

  proc run {steps} {
    variable atoms
    variable groups
    variable script
    variable results
    variable step
    variable oname
    variable f1
    variable f2

    # Run tclforcesscript 
    namespace eval :: {
      source $mock::script
    }

    for {set f 0} {$f < $steps} {incr f} {
        calcforces

        foreach g $groups {
            puts $f1 $results($g.force)
        }            
        puts $f2 $results(energy)

        incr step 
    }
  }
        
  # Auxiliary functions

  proc loadxyz {filename} {
      variable atoms

      set f [open $filename r]

      set nf 1
      while {![eof $f]} {
          if {[gets $f atoms(na)] < 0} { break }
          if {[gets $f comment] < 0} { break }

          for {set i 1} {$i <= $atoms(na)} {incr i} {
              gets $f xyzline
              set fields [regexp -inline -all {\S+} $xyzline]
              set atoms(a$i.sym) [lindex $fields 0]
              set atoms(f$nf.a$i.pos) [lrange $fields 1 4]
          }
          incr nf
      }

      set atoms(nf) $nf
      close $f
  }
       
  # exporta todo (menos las variables)
  namespace export *
         
}

 
namespace import ::mock::*

mock::make_stubs [lindex $argv 0]

source "[lindex $argv 0]"
         
