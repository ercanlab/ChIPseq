#!/bin/bash

for f in reports/*; do
  error="$(cat $f | grep 'error')"
  if [[ -n $error ]]; then
   printf "## File $f has the following errors:\n\t$error\n"
   echo "----------"
   echo
  fi
done
