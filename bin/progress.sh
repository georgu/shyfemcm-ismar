#!/bin/bash
#
# progress bar implementation
#
# Usage: ProgressBar currentState totalState
#
#------------------------------------------------------------------

function ProgressBar {

    # Process data

    let _progress=(${1}*100/${2}*100)/100
    let _done=(${_progress}*4)/10
    let _left=40-$_done

    # Build progressbar string lengths

    _fill=$(printf "%${_done}s")
    _empty=$(printf "%${_left}s")

    # 1.2 Build progressbar strings and print the ProgressBar line
    # 1.2.1 Output example:                           
    # 1.2.1.1 Progress : [########################################] 100%

    printf "\rProgress : [${_fill// /\#}${_empty// /-}] ${_progress}%%"
}

#------------------------------------------------------------------

test_progress()
{
  _start=1
  _end=100

  for number in $(seq ${_start} ${_end})
  do
    sleep 0.1
    ProgressBar ${number} ${_end}
  done
  printf '\nFinished!\n'
}

#------------------------------------------------------------------

#-------------------------------------------------
(return 0 2>/dev/null) || test_progress
#-------------------------------------------------

