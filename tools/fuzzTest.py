#!/usr/bin/env python3
import random, math, functools


def generateExample(nofDimensions,nofPoints,seed,outFile):

    random.seed(seed)
    
    # Generate points on sphere
    points = []
    for i in range(0,nofPoints):
        point = [random.gauss(0,1) for i in range(nofDimensions)]
        pointDist = math.sqrt(sum([a*a for a in point]))
        point = [a/pointDist for a in point]
        points.append(point)

    # Rotate and Scale
    volume = nofDimensions/2.0*math.pow(math.pi,nofDimensions/2.0)/math.gamma(nofDimensions/2.0);
    for i in range(0,50):
        # Scale
        scaleDim = random.sample(range(nofDimensions),1)[0]
        scale = random.gauss(1,0.3)
        for j,point in enumerate(points):
            point[scaleDim] *= scale
        volume *= scale
        
   
        # Rotate
        rotDim1,rotDim2 = random.sample(range(nofDimensions),2)
        angle = random.random()*2*3.1415
        cosAngle = math.cos(angle)
        sinAngle = math.sin(angle)
        for j,point in enumerate(points):
            m1 = point[rotDim1]
            m2 = point[rotDim2]
            c1 = m1*cosAngle - m2*sinAngle
            c2 = m1*sinAngle + m2*cosAngle
            point[rotDim1] = c1
            point[rotDim2] = c2

    with open(outFile,"w") as outf:
        for p in points:
            outf.write(" ".join([str(a) for a in p]))
            outf.write("\n")

    return volume


volume = generateExample(10,1000,10,"/tmp/test.txt")
print("Volume: ",volume)

