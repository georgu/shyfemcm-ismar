#!/bin/bash
#
# progress bar implementation
#
# Usage: ProgressBar currentState totalState
#
#------------------------------------------------------------------

stringProgressBar="Progress: "

function ProgressBar {

    # Process data

    [ "$2" = "0" ] && return
    [ "$1" -gt "$2" ] && return

    let _progress=(${1}*100/${2}*100)/100
    let _done=(${_progress}*4)/10
    let _left=40-$_done
    local _string="$stringProgressBar"

    # Build progressbar string lengths

    _fill=$(printf "%${_done}s")
    _empty=$(printf "%${_left}s")

    # 1.2 Build progressbar strings and print the ProgressBar line
    # 1.2.1 Output example:                           
    # 1.2.1.1 Progress : [########################################] 100%

    printf "\r${_string} [${_fill// /\#}${_empty// /-}] ${_progress}%%"
}

function initProgressBar {

    stringProgressBar=$1
}

function finalizeProgressBar_with_nl {

  finalizeProgressBar
  echo ""
}

function finalizeProgressBar {

    let _progress=100
    let _done=40
    let _left=0
    local _string="$stringProgressBar"

    _fill=$(printf "%${_done}s")
    _empty=$(printf "%${_left}s")

    printf "\r${_string} [${_fill// /\#}${_empty// /-}] ${_progress}%%"
}

#------------------------------------------------------------------

test_progress()
{
  _start=1
  _end=100

  initProgressBar "my_string"

  for number in $(seq ${_start} ${_end})
  do
    sleep 0.1
    ProgressBar ${number} ${_end}
  done
  finalizeProgressBar
  printf ' Finished!\n'
}

#------------------------------------------------------------------

#-------------------------------------------------
(return 0 2>/dev/null) || test_progress
#-------------------------------------------------

