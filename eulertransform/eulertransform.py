##
## Comparison of the linear and the quadratic implementation of the Euler
## Tranform, applied to finding decimals of pi.
## For details on the algorithm, see 
## https://rmj86.github.io/blog/2016/11/19/euler-transform-in-linear-time
##

from decimal import *

def pi_series(n):
    """ return the first n terms in leibnit's seies for pi
        (alternating sign) """
    elems = [(-1)**i * Decimal(4)/Decimal(2*i+1) for i in xrange(n)]
    return elems

def _nabla(series, n):
    """ nth term in euler transform of series. Use to estimate precision. 
        Provided that the terms decrease quick enough, the error in
        partialT(series, n) will be no larger than _nabla(series, n)."""
    nab = Decimal(0)
    for a, f in zip(series, pasc_row_lin(n)):
        nab += a * f
    return nab / (2**(n+1))

################################################################################
## Standard formula    O(n^2)

def next_pasc_row(r):
    next_row = [1]
    for a,b in zip(r, r[1:]):
        next_row.append(a+b)
    next_row.append(1)
    return next_row

def partialT(series, n):
    terms = []
    pasc_row = [1]
    for k in xrange(0,n+1):
        nab = Decimal(0)
        for a, f in zip(series, pasc_row):
            nab += a * f
        terms.append(nab / (2**(k+1)))
        pasc_row = next_pasc_row(pasc_row)
    s = sum(reversed(terms))
    return s
    

################################################################################
## simplified formula    O(n)

def pasc_row_lin(n):
    b = Decimal(1)
    row = [b]
    for i in range(n):
        b *= Decimal(n-i) / Decimal(i+1)
        row.append(b)
    return row

def et_consts(n):
    """ constant terms for the nth partial sum of the simplified ET formula """
    c = Decimal(2**(n+1))
    consts = []
    for b in pasc_row_lin(n+1):
        c -= b
        consts.append(c)
    return consts[:-1]

def partialT_lin(series, n):
    sum = Decimal(0)
    for a,c in reversed(zip(series, et_consts(n))):
        sum += a*c
    return sum / Decimal(2**(n+1))

################################################################################

if __name__ == "__main__":
    from time import time
    
    # We get roughly 105 digits of precision for every
    # 350 terms of the series. (empirical for  n_terms <= 1000)
    n_terms = 1 * 350
    getcontext().prec = 1 * 110  # set precision of decimal arithmetic
    
    series = pi_series(n_terms+2)

    # print _nabla(series, n_terms)  ## expected error, roughly

    print "Calculating pi via Euler Transform on leibniz's formula,"
    print "pi = 1 - 1/3 + 1/5 - 1/7 + ..."

    t0 = time()
    print "Using linear formula:"
    print partialT_lin(series, n_terms)
    print "time: {}".format(time() - t0)

    t0 = time()
    print "Using quadratic formula:"
    print partialT(series, n_terms)
    print "time: {}".format(time() - t0)
    print
