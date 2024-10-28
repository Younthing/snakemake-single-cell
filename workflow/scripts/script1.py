# Python script to print "Hello Python"
import os

# Ensure the directory exists
os.makedirs("results", exist_ok=True)

with open("results/hello_python.txt", "w") as f:
    f.write("Hello Python\n")
