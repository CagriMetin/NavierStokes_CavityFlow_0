import numpy
from matplotlib import pyplot
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

x,y,p = numpy.loadtxt('pressure.txt', delimiter=',', unpack=True)
x,y,u = numpy.loadtxt('u_vel.txt', delimiter=',', unpack=True)
x,y,v = numpy.loadtxt('v_vel.txt', delimiter=',', unpack=True)
n = int(numpy.sqrt(len(x)))
assert n*n == len(x), '"Expected len(x) to be a perfect square, len(x) = %s" % len(x) '

X = x.reshape(n,n)
Y = y.reshape(n,n)
P = p.reshape(n,n)
U = u.reshape(n,n)
V = v.reshape(n,n)

fig = pyplot.figure(figsize=(11,7), dpi=100)
pyplot.contourf(X,Y,P,alpha=0.5)    ###plnttong the pressure field as a contour
pyplot.colorbar()
pyplot.contour(X,Y,P)               ###plotting the pressure field outlines
pyplot.quiver(X[::2,::2],Y[::2,::2],U[::2,::2],V[::2,::2]) ##plotting velocity
pyplot.xlabel('X')
pyplot.ylabel('Y');
pyplot.show()
