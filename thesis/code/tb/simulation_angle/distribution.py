import numpy as np
import math

"""
   draw from a 1D distribution using the inverted cdf approach.
   Partionally inspired by: http://stackoverflow.com/questions/21100716/fast-arbitrary-distribution-random-sampling

   The mapping from fractional indices to xvalues is achieved by
   a static double lambda function, very nice! ((c) Camille Colle 2015)
"""

class Distribution(object):
    # xmap maps the fractional indices to original xvalues, supply x array as xr
    xmap = staticmethod(lambda xr: lambda idx: [ (1.-f)*xr[int(i)]+f*xr[int(i)+1] for f,i in map(math.modf,idx) ])

    def __init__(self,pdf,interpolation=True,transform=lambda x:x):
        assert(np.all(pdf>=0))
        self.pdf = pdf
        self.cdf = np.cumsum(self.pdf)
        self.interpolation = interpolation
        self.transform = transform
    def sum(self):
        return self.cdf[-1]
    def __call__(self,N):
        """  draw
             this method returns fractional indices
             you can map this using your own transform
             function.
             If you want to map the indices to your original
             xrange initialize the object with
             dist(self,pdf,interpolation,transform=Distribution.xmap(xr))
             where xr is the xrange you are using.
             This will likely be the most used option
        """
        choice = np.random.uniform(high=self.sum(),size=N)
        index  = np.searchsorted(self.cdf,choice)
        if self.interpolation:
            index = index+np.random.uniform(size=N)
        return self.transform(index)

if __name__ == "__main__":
    import scipy.stats
    xr  = np.linspace(-5.,5.,1001)
    pdf = scipy.stats.norm.pdf(xr) 
    dist = Distribution(pdf,transform=Distribution.xmap(xr))

    import pylab as pl
    fig = pl.figure()
    ax = fig.add_subplot(111)
    ax.hist(dist(100000),normed=True,bins=21)
    ax.plot(xr,pdf,'r--',lw=3)
    pl.show()


