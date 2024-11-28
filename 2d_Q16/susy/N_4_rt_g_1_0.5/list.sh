#!/bin/bash

# Define the maximum value for the range
n=2000  # You can change this to your desired value

# Initialize the output file
output_file="list.txt"
> "$output_file"

# Generate the ranges and write to the file
for ((i=0; i < n; i+=10)); do
    start=$i
    end=$((i+10))
    echo "${start}-${end}" >> "$output_file"
done

# Add the final range if necessary
if ((n % 10 != 0 || n >= 10)); then
    echo "${((n/10)*10)}-${n}" >> "$output_file"
fi

echo "list.txt has been created with ranges from 0-10 to $(($n/10*10))-${n}."

