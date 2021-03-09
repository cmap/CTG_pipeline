#!/bin/bash

# read in flagged arguments
while getopts ":f:m:r:o:p:" arg; do
  case $arg in
    f) # specify mapping path
      merged=${OPTARG};;
    m) # specify mergedfile path
      mapping=${OPTARG};;
    r) # specify raw data path
      raw=${OPTARG};;
    o) # specifcy output folder
      output_dir=${OPTARG};;
    p) # specify the name of the project
      project=${OPTARG};;
  esac
done

Rscript ./process_CTG.R "${merged}" "${mapping}" "${raw}" "${output_dir}" "${project}"

exit_code=$?

echo "$exit_code"
exit $exit_code
