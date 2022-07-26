#!/usr/bin/env python3
import os, sys, numpy

if len(sys.argv)<3:
    print("Error: Expected two parameters:\n- Output File of Elipsoid Computation Tool.\n- Test Vector File\n")
    sys.exit(1)
    
centerVector = None
matrix = []

# Read input file with Elipsoid    
with open(sys.argv[1],"r") as inFile:
    nextCenterVector = False
    nextMatrix = False
    for line in inFile.readlines():
        line = line.strip()
        if line=="Result - Center vector:":
            nextCenterVector = True
        elif line=="Result - Matrix:":
            nextMatrix = True
        elif line.startswith("Result"):
            nextMatrix = False
        elif nextCenterVector:
            centerVector = [float(a) for a in line.split(" ")]
            nextCenterVector = False
        elif nextMatrix:
            matrix.append([float(a) for a in line.split(" ")])

# Sanity check that the dimensions match
print("# Dimensions Elipsoid: ",len(centerVector))
if len(matrix)!=len(centerVector):
    print("Error: Dimension of matrix does not match (A)")
    sys.exit(1)
for a in matrix:
    if len(a)!=len(centerVector):
        print("Error: Dimension of matrix does not match (B)")
        sys.exit(1)


# Translate to numpy
nofDims = len(centerVector)
matrix = numpy.array(matrix)
centerVector = numpy.array(centerVector)

# Try the data points
with open(sys.argv[2],"r") as inFile:
    for line in inFile.readlines():
        line = line.strip()
        if len(line)>0:
            thisVector = [float(a) for a in line.split(" ")]
            if len(thisVector)<len(centerVector):
                print("Error: Dimension mismatch in points file")
                sys.exit(1)
            thisVector = numpy.array(thisVector[0:nofDims])
            
            # Compute distance
            dist = thisVector - centerVector
            # print(dist.shape)
            distT = numpy.transpose([dist])
            # print(distT.shape)
            dist = numpy.matmul(dist,numpy.matmul(matrix,distT))
            print(dist[0])
            

