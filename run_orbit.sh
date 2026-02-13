#!/bin/bash

# Step 1: Compile the program
make

# Step 2: Run the compiled executable
./orbit

# Step 3: Clean up the build files
make clean

# Step 4: Run the Python post-processing script
# python eo_post_process_cpp.py

# Step 5 delete the output file
# rm output.csv