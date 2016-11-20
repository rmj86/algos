##
## Comparison of the linear and the quadratic implementation of the Euler
## Tranform, applied to finding decimals of pi.
## For details on the algorithm, see
## https://rmj86.github.io/blog/2016/11/19/euler-transform-in-linear-time
##

from decimal import *

def leibniz_pi(n):
    """ return the first n terms in Leibniz's seies for pi
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
    for i in range(len(r)-1):
        next_row.append( r[i]+r[i+1] )
    next_row.append(1)
    return next_row

def partialT(series, n):
    """ nth partial sum of series' Euler Transform. O(n^2) """
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
    """ Row n of pascal's triangle. O(n) """
    b = 1
    row = [b]
    for i in range(n):
        b = (b*(n-i)) / (i+1)
        row.append(b)
    return row

def et_consts(n):
    """ constant terms for the nth partial sum of the simplified ET formula """
    c = 2**(n+1)
    consts = []
    for b in pasc_row_lin(n+1):
        c -= b
        consts.append(c)
    return consts[:-1]

def partialT_lin(series, n):
    """ nth partial sum of series' Euler Transform. O(n) """
    sum = Decimal(0)
    for a,c in reversed(zip(series, et_consts(n))):
        sum += a*c
    return sum / Decimal(2**(n+1))


################################################################################

if __name__ == "__main__":
    from time import time

    # We get roughly 105 digits of precision for every
    # 350 terms of the series. (empirical for  n_terms <= 1000)
    digits = 100
    n_terms = int(digits*3.5)
    getcontext().prec = digits+10  # set precision of decimal arithmetic

    series = leibniz_pi(n_terms+2)


    print "Calculating pi via Euler Transform on leibniz's formula,"
    print "pi = 1 - 1/3 + 1/5 - 1/7 + ... \n"
    # print "{} terms, expected error {:.5e}\n".format(n_terms, 10*_nabla(series, n_terms))
    
    t0 = time()
    print "Using linear formula:"
    print partialT_lin(series, n_terms)
    print "time: {}\n".format(time() - t0)
    
    t0 = time()
    print "Using quadratic formula:"
    print partialT(series, n_terms)
    print "time: {}\n".format(time() - t0)
    