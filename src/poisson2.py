#!/usr/bin/env python
import sys
sys.path.append("/home/ataias/workspace/ff/build/src/")

from libfooclass import FooClass
import numpy
m = 7
n = 8
xIn = numpy.ones([m,n])
xIn[1][1]=0
xIn[2][1]=2.954
xOut = numpy.zeros([m,n])

print "Before"
print xIn
print xOut

f = FooClass(m,n)

f.foo(xIn, xOut)

print xIn
print xOut