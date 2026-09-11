import numpy as np

xlow = np.log10(0.002)
xhigh = np.log10(np.pi/2)
nbins = 100

width = (xhigh-xlow)/nbins

bins=[]
for i in range(nbins+1):
    val = pow(10, xlow + i * width)
    bins += [val]

newbins = [np.pi- b for b in bins]
newbins = newbins[::-1]
del newbins[0]

bin_edge = bins+newbins

print(bin_edge)
