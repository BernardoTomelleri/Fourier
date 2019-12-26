# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 23:03:46 2019

@author: berna
"""

import numpy
import scipy.integrate as integrate

def scalar_product(f, g, domain):
    ans = integrate.quad( lambda x: f(x)*g(x), *domain )[0]
    return ans

def norm(f, domain):
    return scalar_product(f, f, domain)
    
"""
Generico spazio L2. Parametri:
-domain: dominio dello spazio
"""    
class L2Space:
    def __init__(self, domain):
        self.domain = domain

class FourierBase:
    def __init__(self, L2, f):
        self.L2 = L2
        self.f = f
        
    def e(self, index, x, offset = None):
        return self.f(x, index, self.L2.domain, offset)
    
    def coefficient(self, index, function):
        return scalar_product( function, lambda x: self.f(x, index, self.L2.domain), self.L2.domain )
        
    def coefficents(self, f, index):
        i = index.zero()
        ans = []
        while i <= index:
            ans += [self.coefficient(i, f)]
            i.increment()
        return ans
    
    def fcoefficients(self, coeff, i0, x):
        ans = numpy.zeros(len(x))
        for k in range(0,len(coeff)):
            #print(i0, coeff[k])
            ans += coeff[k] * FourierBase.e(self,i0, x)
            i0.increment()
        return ans
    
    def fseries(self, f, index, x):
        return self.fcoefficients(self.serie(f, index), x, index.zero())
    
    

class cossinindex:
    def __init__(self, n=0, f = 'cos'):
        self.n = n
        if f != 'cos' and f != 'sin':
            raise TypeError('invalid type {}'.format(f))
        self.f = f
        
    def increment(self):
        if self.n == 0:
            self.n += 1
        elif self.f == 'cos':
            self.f = 'sin'
        elif self.f == 'sin':
            self.f = 'cos'
            self.n += 1
        else:
            raise IndexError('invalid index ({0},{1})'.format(self.n, self.f))
            
    def zero(self):
        return cossinindex()
        
    def __eq__(self, other):
        return self.n == other.n and self.f == other.f
        
    def __lt__(self, other):
        return self.n < other.n
        
    def __le__(self, other):
        return self.n <= other.n
    
    def __repr__(self):
        return '({0}, {1})'.format(self.n, self.f)
       
def cossinnormcoefficeint(index, domain):
    numpy.sqrt( 2 / (domain[1]-domain[0]))

def cossinbasefunction(x, index, domain, offset = None):
    if offset == None:
        offset = numpy.zeros(len(x))
    if index.f == 'sin':
        if index.n != 0:
            return numpy.sin(2*numpy.pi*index.n/(domain[1]-domain[0]) * (x-domain[0]) + offset) * numpy.sqrt( 2 / (domain[1]-domain[0]))
        else:
            raise IndexError('invalid index for sincosbase {}'.format(index))
    elif index.f == 'cos':
        if index.n == 0:
            return numpy.cos(0.) * numpy.sqrt(1/ (domain[1]-domain[0]))
        else:
            return numpy.cos(2*numpy.pi*index.n/(domain[1]-domain[0]) * (x-domain[0]) + offset) * numpy.sqrt( 2 / (domain[1]-domain[0]))
    

class CosSinBase(FourierBase):
    def __init__(self, L2, os = False, oc = False):
        super().__init__(L2, cossinbasefunction)
        
    def e(self, order, x, offset = None):
        if order != 0:
            return [super().e(cossinindex(order, 'cos'), x, offset), super().e(cossinindex(order, 'sin'), x, offset)]
        return [super().e(cossinindex(order, 'cos'), x, offset), None]
        
    def coefficents(f, order):
        ans = [super(FourierBase, self).coefficient(cossinindex(), f), None]
        for i in range(1, order+1):
            ans += [super().coefficient(cossinindex(i, 'cos'), f), super().coefficient(cossinindex(i, 'sin'), f)]
        return ans
        
    def fseries(self,f, order, x):
        return super().fseries(f, cossinindex(order), x)
        
    def fcoefficients(self, coeff, x):
        return super().fcoefficients( [coeff[0][0]] + [ coeff[c][it] for c in range(1,len(coeff)) for it in range(0,2) ], cossinindex(), x )


class SinBase(CosSinBase):
    def e(self, order, x, offset = None):
        if order != 0:
            return super().e(order,x, offset)[1]
        return super().e(order,x, offset)[0]
    
    def coefficents(self, f, order):
        sans = super().coefficents(f, order)
        return [sans[0][0]] + [sans[i][1] for i in range(1,order+1)]
        
    def fcoefficients(self, coeff, x):
        return super().fcoefficients([[coeff[0], None]] + [ [0., coeff[i]] for i in range(1, len(coeff))], x)
    
class CosBase(CosSinBase):
    def e(self, order, x, offset = None):
        return super().e(order, x, offset)[0]
        
    def coefficents(self, f, order):
        sans = super().coefficents(f, order)
        return [sans[0][0]] + [sans[i][0] for i in range(1,order+1)]
        
    def fcoefficients(self, coeff, x):
        return super().fcoefficients([[coeff[0], None]] + [ [coeff[i], 0.] for i in range(1, len(coeff))], x)
    
