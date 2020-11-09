#!/bin/csh -f

# Don't source .cshrc file.  Modules may be swapped/loaded confusing the batch scheduler

#______________________________________________________________________
#  monitorSus
#  This script searches the output file every N minutes and 
#  if the timestep number has not increased then the job ID is killed.
#______________________________________________________________________

getopt -Q -q --long pid,out,interval -- $argv:q >& /dev/null

if($#argv != 6) then
  echo ""
  echo "Usage: $0 "
  echo "This script searches the output file every N minutes for"
  echo " an increase in the timestep number.  If the timestep"
  echo " is not increasing the jobID is killed with a scancel command"
  echo ""
  echo " Manditory inputs:"
  echo "   --jid      <job id that will be killed>"
  echo "   --out      <name of the outputfile to monitor>"
  echo "   --interval <interval in minutes between searches>"
  exit 1
endif

while ($#argv)
   switch ($1:q)
     case --jid:
        set jid = $2
        shift; shift;
        breaksw
     case --out:
        set out = "$2"
        shift; shift
        breaksw
     case --interval:
        set interval = "$2"
        shift; shift
        breaksw
     case " ":
        shift
        break
     default:
        echo "Usage: $0 "
        echo ""
        echo " Manditory inputs:"
        echo "   --jid      <job id that will be killed>"
        echo "   --out      <name of the outputfile to monitor>"
        echo "   --interval <interval in minutes between searches>"
        exit 1
   endsw
end

set tmp = (`which scancel`)
if ( $status ) then
  echo ""
  echo " ERROR: Could not find the command scancel."
  echo ""
  exit(-1)
endif


@ ts_old = 0

while (1)
  set minutes = $interval"m"
  sleep $minutes
  
  @ ts = `grep Timestep $out | tr -d '[:punct:]' | awk 'END {print $2}'`
  if $ts == $ts_old then
    date
    echo "  MonitorSus.csh: Now killing job $jid on timestep $ts_old"
    scancel $jid
    
    exit(1)
  else
    @ ts_old = $ts
  endif
end

exit 
