#!/usr/bin/tclsh

topo

topo clearbonds

set sel [atomselect top "name 1 or name 2"]

set Ntotal [$sel num ]
 
set Nmon 100

puts $Ntotal

set Ntotal [expr $Ntotal/$Nmon]

puts $Ntotal

for {set a 0} {$a < $Ntotal} {incr a} {
   for { set i [expr $a*$Nmon]} {$i < [expr ($a+1)*$Nmon-1]} {incr i} { 
      set j [expr $i+1]	
      topo addbond $i $j 
   }
} 



