# Python script to print "Hello Python"
import os
import sys

# Ensure the directory exists
os.makedirs("results", exist_ok=True)

output_file = sys.argv[1]

with open(output_file, "w") as f:
    f.write("Hello Python\n")
