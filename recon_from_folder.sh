#!/bin/bash
args_array=("$@")
for i in "${args_array[@]}"
do
  :
  echo "variable $i ###"
done
echo "args_count = $#"
# join the arguments as single string to execute as command
eval "$@"	
