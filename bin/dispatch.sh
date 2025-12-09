#!/bin/bash
#
# dispatch routines for running jobs in batch
#
#---------------------------------------------------

debug="YES"
debug="NO"
date=$( date +%Y-%m-%d::%H:%M:%S )
dateaux=$( date +%Y-%m-%d__%H_%M_%S )
nstats=10
nmin=$(( nstats / 2 ))
hpdir=~/tmp/hp
hpbatch=$hpdir/hpbatch
hplog=$hpdir/hplog
hplogfile=$hplog/batch.$dateaux.log
logfile=cpu.log
secs=2

mkdir -p $hplog
mkdir -p $hpbatch
mkdir -p $hpbatch/finished

#---------------------------------------------------

RunCpuStats()
{
  [ $# -eq 0 ] && exit 1
  local nwanted=$1

  while true
  do
    sleep $secs
    GetCpuStats $nwanted
    date=$( date +%Y-%m-%d::%H:%M:%S )
    #echo "$date $cpufree $npositive $positive" | tee -a $hplogfile
    echo "$date $cpufree $npositive $positive" >> $hplogfile
    [ $positive = "YES" ] && break
  done
}

GetCpuStats()
{
  [ $# -eq 0 ] && return
  local nwanted=$1

  npositive=0

  for i in $( seq $nstats )
  do
    GetCpuUsage $nwanted
    [ $canrun = "YES" ] && npositive=$(( npositive + 1 ))
  done

  [ $debug = "YES" ] && echo "positive: $npositive $nmin"

  positive="NO"
  [ $npositive -ge $nmin ] && positive="YES"
}

GetCpuUsage()
{

canrun="NO"
[ $# -eq 0 ] && return
local nwanted=$1

nproc=$( nproc )
pwanted=$(( nwanted * 100 / nproc ))

line=$( top -b -n 1 | grep Cpu | sed -e 's/:/ /' | sed -e 's/,/ /g' )
[ $debug = "YES" ] && echo $line

totcpu=0
for i in 2 4 6
do
  cpu=$( echo $line | cut -d " " -f $i | sed -e 's/\..*//' )
  totcpu=$(( totcpu + cpu ))
done

cpufree=$(( 100 - totcpu ))

if [ $debug = "YES" ]; then
  echo "total cpu: $totcpu"
  echo "total number of processors: $nproc"
  echo "processors wanted: $nwanted"
  echo "percentage wanted: $pwanted"
  echo "free cpu: $cpufree"

  if [ $cpufree -ge $pwanted ]; then
    echo "can run new process"
  else
    echo "cpus are occupied... cannot run"
  fi
fi

[ $cpufree -ge $pwanted ] && canrun="YES"

[ $debug = "YES" ] && echo "$nproc $nwanted $totcpu $cpufree $pwanted $canrun"

}

#---------------------------------------------------
#---------------------------------------------------
#---------------------------------------------------

#Dispatch_loop_until_run()
#{
#}

Dispatch_run_batch()
{
  local file=$1
  echo "handling file $file"
  np=$( grep 'npmpi=' $file | sed -e 's/npmpi=//' )
  [ $np -eq 0 ] && np=1
  echo "np = $np"
  RunCpuStats $np
  #Dispatch_loop_until_run $np
  #exit 0
  echo "running file $file"
  bash $file		!run simulation
  sleep 10
}

Dispatch_is_running()
{
  local nn=$( lsof -t $0| wc -l )

  echo "running processes $nn"

  if [ $nn -gt 1 ]; then
    return 0
  else
    return 1
  fi
}

Dispatch_start()
{
  if Dispatch_is_running; then
    echo "dispatch is already running... exiting"
    return
  fi

  echo "starting dispatch"
  echo "starting dispatch" >> $hplogfile

  local n=0
  while true
  do
    sleep 4
    file=$( ls -1 $hpbatch/hp* 2> /dev/null | head -1 )
    #echo $file
    if [ -z "$file" ]; then
      echo "no file left... exiting" >> $hplogfile
      break
    fi
    echo "dealing with $file" >> $hplogfile
    # here we have to run the command when possible
    Dispatch_run_batch $file
    mv -f $file $hpbatch/finished
    n=$(( n + 1 ))
    [ $n -ge 6 ] && break
  done
}

Dispatch_stop()
{
  echo "stopping dispatch" >> $hplogfile
  exit 0
}

Dispatch_test()
{
  echo "starting dispatch test"

  #touch $hpbatch/hp.1 $hpbatch/hp.2 $hpbatch/hp.3

  Dispatch_start

  echo "finished dispatch test"

  Dispatch_stop

  if Dispatch_is_running; then
    echo "dispatch is still running"
  fi
}

#---------------------------------------------------

#Dispatch_test
Dispatch_start

#---------------------------------------------------

