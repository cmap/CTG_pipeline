#!/bin/bash

# read in flagged arguments
while getopts ":i:o:p:" arg; do
  case $arg in
    i) # specify input folder
      data_dir=${OPTARG};;
    o) # specifcy output folder
      output_dir=${OPTARG};;
    p) # specify the name of the project
      project=${OPTARG};;
  esac
done

Rscript /process_CTG.R "${data_dir}" "${output_dir}" "${project}"

exit_code=$?

echo "$exit_code"
exit $exit_code
