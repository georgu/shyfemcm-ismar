#!/bin/bash
#
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------
#
# examples and usage of advanced features of bash
#
#------------------------------------------------------------------------

# uncomment next command to catch errors in script
# notice: some commands will not run

#set -euo pipefail

#------------------------------------------------------------------------

Test_Var()
{
  # https://www.howtogeek.com/782514/

  set -e 		# exit after error
  set -eo pipefail	# exit after error, even in pipes
  set -u 		# complain because of non initialized variables
  set -x		# write execution of script
  set -euxo pipefail	# all of the above
}

Test_Func()
{
  return $1
}

Test_If()
{
  echo "================================="
  echo "test if"
  echo "================================="

  if Test_Func 0; then echo "ok"; else echo "false"; fi
  if Test_Func 1; then echo "ok"; else echo "false"; fi
}

Test_PatternMatching()
{
  echo "================================="
  echo "pattern matching"
  echo "================================="

  vech="Bar"
  echo "match for B* in $vech:"
  [[ $vech = B* ]] && echo "Start with B" || echo "Not matched"
  vech="Car"
  echo "match for B* in $vech:"
  [[ $vech = B* ]] && echo "Start with B" || echo "Not matched"
}

Test_Loop()
{
  echo "================================="
  echo "test loop"
  echo "================================="

  echo "----------------Loop over range---------------------"
  for id in {1..5};
  do
    echo "$id"
  done
}

Test_ReadFromFile()
{
  echo "================================="
  echo "read from file"
  echo "================================="

  local file=$0
  echo "reading from file $file"

  local ilines=0
  while read -r line; do
    (( ilines++ ))
    #echo "$line"
  done < $file

  echo "$ilines read"
}

Test_Arithmetics()
{
  echo "================================="
  echo "arithmetics"
  echo "================================="

  a=$(( 4 + 5 ))
  echo "$a == 9"
  a=$((3+5))
  echo "$a == 8"
  b=$(( a + 3 ))
  echo "$b == 11"
  b=$(( $a + 4 ))
  echo "$b == 12"
  (( b++ ))
  echo "$b == 13"
  (( b += 3 ))
  echo "$b == 16"
  a=$(( 4 * 5 ))
  echo "$a == 20"
}

Test_Array()
{
  echo "================================="
  echo "test array"
  echo "================================="

  MyArray=( HTML Javascript CSS JQuery Bootstrap )
 
  echo "----------Print 5 values individually---------------"
  echo ${MyArray[0]}
  echo ${MyArray[1]}
  echo ${MyArray[2]}
  echo ${MyArray[3]}
  echo ${MyArray[4]}
 
  echo "-----------------Print all values-------------------"
  echo ${MyArray[*]}
  echo "-----------------Print all values-------------------"
  echo ${MyArray[@]}

  echo "----------------Loop through values-----------------"
  for i in ${!MyArray[@]}; do
    echo "$i - ${MyArray[$i]}"
  done

  echo "--------------------Add to array--------------------"
  allThreads=(1 2 4 8 16 32 64 128)
  allRuntimes=()
  for t in ${allThreads[@]}; do
    runtime=$t
    allRuntimes+=( $runtime )
  done
  echo ${allThreads[*]}
  echo ${allRuntimes[*]}
  echo "array size: ${#allThreads[@]}"

  echo "------------------Shift first item------------------"
  allThreads=("${allThreads[@]:1}")
  echo ${allThreads[*]}
  echo "array size: ${#allThreads[@]}"

  echo "-------------------Pop last item--------------------"
  allThreads=("${allThreads[@]:0:${#allThreads[@]}-1}")
  echo ${allThreads[*]}
  echo "array size: ${#allThreads[@]}"

  echo "------------------Get items 2,3,4-------------------"
  allThreads=("${allThreads[@]:2:3}")
  echo ${allThreads[*]}
  echo "array size: ${#allThreads[@]}"

  # arr=(1 2 3) 	# Initialize array
  # ${arr[2]} 		# Retrieve third element
  # ${arr[@]} 		# Retrieve all elements
  # ${!arr[@]} 		# Retrieve array indices
  # ${#arr[@]} 		# Calculate array size
  # arr[0]=3 		# Overwrite 1st element
  # arr+=(4) 		# Append value(s)
  # str=$(ls) 		# Save ls output as a string
  # arr=( $(ls) ) 	# Save ls output as an array of files
  # ${arr[@]:s:n} 	# Retrieve n elements starting at index s
}

Test_Associative()
{
  echo "================================="
  echo "test associative array"
  echo "================================="

  declare -A assArray1

  assArray1[fruit]=Mango
  assArray1[bird]=Cockatail
  assArray1[flower]=Rose
  assArray1[animal]=Tiger

  echo "------------------show keys of array-------------------"
  echo "${!assArray1[@]}"
  for key in "${!assArray1[@]}"
  do 
    echo $key
  done

  echo "------------------show vals of array-------------------"
  echo "${assArray1[@]}"
  for val in "${assArray1[@]}"
  do 
    echo $val
  done

  echo "------------------show keys and vals of array-------------------"
  for key in "${!assArray1[@]}"
  do 
    echo "$key: ${assArray1[$key]}"
  done

  echo "------------------access elements-------------------"
  echo "bird: ${assArray1[bird]}"
  echo "flower: ${assArray1[flower]}"

  echo "------------------add new element-------------------"
  echo "${assArray1[@]}"
  assArray1+=([vegatable]=Potato)
  echo "${assArray1[@]}"

  echo "------------------delete element-------------------"
  echo "${assArray1[@]}"
  unset assArray1[fruit]
  echo "${assArray1[@]}"
  echo "fruit: ${assArray1[fruit]}"	# this throws an error

  echo "------------------look for element-------------------"
  echo -n "look for flower: "
  if [ ${assArray1[flower]+_} ]; then echo "Found"; else echo "Not found"; fi
  echo -n "look for fruit: "
  if [ ${assArray1[fruit]+_} ]; then echo "Found"; else echo "Not found"; fi
  echo -n "look for fruit: "	# next command also throws error
  if [[ -v ${assArray1[fruit]} ]]; then echo "Found"; else echo "Not found"; fi

  echo "------------------delete array-------------------"
  echo "${assArray1[@]}"
  unset assArray1
  echo "${assArray1[@]}"
}

Test_Round()
{
  echo "================================="
  echo "test rounding"
  echo "================================="

  val="2.3"; new=$( round $val 0 ); echo "$val -> $new (rounding to pos 0)"
  val="2.5"; new=$( round $val 0 ); echo "$val -> $new (rounding to pos 0)"
  val="2.7"; new=$( round $val 0 ); echo "$val -> $new (rounding to pos 0)"
  val="-2.3"; new=$( round $val 0 ); echo "$val -> $new (rounding to pos 0)"
  val="-2.5"; new=$( round $val 0 ); echo "$val -> $new (rounding to pos 0)"
  val="-2.7"; new=$( round $val 0 ); echo "$val -> $new (rounding to pos 0)"
  val="2.34"; new=$( round $val 1 ); echo "$val -> $new (rounding to pos 1)"
  val="2.56"; new=$( round $val 1 ); echo "$val -> $new (rounding to pos 1)"
  val="2.71"; new=$( round $val 1 ); echo "$val -> $new (rounding to pos 1)"
}

Test_Messages()
{
  echo "================================="
  echo "test messages (should be last test)"
  echo "================================="

  ok This is an ok message
  info This is an info message
  warn This is a warning message
  die This is an error

  info This message will never be shown
}

Test_Directory()
{
  if Exists_directory $1; then ok "directory $1 exists"; \
		else warn "directory $1 does not exist"; fi
}

Test_File()
{
  if Exists_file $1; then ok "file $1 exists"; \
		else warn "file $1 does not exist"; fi
}

Test_Delete()
{
  if Delete_file $1; then ok "file $1 has been deleted"; \
		else warn "file $1 could not be deleted"; fi
}

Test_Files()
{
  echo "================================="
  echo "test files"
  echo "================================="

  Test_Directory tmp
  Test_Directory tmptmp
  Test_File bashutil.sh
  Test_File bashutil.sh7
  touch tmp.tmp
  Test_File tmp.tmp
  Test_Delete tmp.tmp
  Test_Delete tmp.tmp
}

#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------
# utility routines
#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------

ok() { echo -e "\033[1;32m[OK]\033[0m $*"; }
info() { echo -e "\033[1;34m[INFO]\033[0m $*"; }
warn() { echo -e "\033[1;33m[WARN]\033[0m $*"; }
die() { echo -e "\033[1;31m[ERRORE]\033[0m $*" >&2; exit 1; }

#-----------------------------------------------------------

Delete_file()
{
  [ $# -eq 0 ] && return 1
  [ ! -f $1 ] && return 1
  rm -f $1
}

Exists_file()
{
  [ $# -eq 0 ] && return 1
  [ ! -f $1 ] && return 1
  return 0
}

Exists_directory()
{
  [ $# -eq 0 ] && return 1
  [ ! -d $1 ] && return 1
  return 0
}

#-----------------------------------------------------------

log()
{ 
  # returns floor[log(x)] for x integer
  # call as "log 777" (base 2) or "log 10 777" (base 10)
  # http://phodd.net/gnu-bc/bcfaq.html#bashlog

  local x=$1 n=2 l=-1;
  if [ "$2" != "" ]; then 
    n=$x;x=$2
  fi

  while((x))
  do
    let l+=1 x/=n
  done

  echo $l;
}

round()
{
  # call as "round var [ipos]" : rounds var to ipos positions (default ipos=2)

  local val=$1
  local ipos=$2
  [ -z "$ipos" ] && ipos=2

  printf "%.${ipos}f" "${val}"
}

#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------
# execute Test routines
#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------

Test()
{
  Test_Loop
  Test_Array
  Test_Associative
  Test_If
  Test_Arithmetics
  Test_ReadFromFile
  Test_PatternMatching
  Test_Round
  Test_Files
  Test_Messages
}

#-----------------------------------------------------------

(return 0 2>/dev/null) || Test

#-----------------------------------------------------------

