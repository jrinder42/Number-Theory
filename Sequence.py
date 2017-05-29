import numpy as np
import math
import mpmath # nsum, inf
import scipy.misc
from Primes import iscoprime, prime_factorization, divisors

def karatsuba(x,y): # FIX so it is not recursion
    pass

# (Private) Returns the tuple (F(n), F(n+1)).
def _fib(n):
    # Helper function for fibonacci_n(n)
    if n == 0:
        return (0, 1)
    else:
        a, b = _fib(n // 2)
        c = a * (b * 2 - a)
        d = a * a + b * b
        if n % 2 == 0:
            return (c, d)
        else:
            return (d, c + d)

def fibonacci_n(n): # O(log n)
    '''
    :param n: nth fibonacci number requested
    :return: nth fibonacci number
    '''
    if n < 0:
        raise ValueError("Negative arguments not implemented")
    return _fib(n)[0]


def fibonacci(n):

    return (fibonacci_n(i) for i in xrange(n))

def lucas_n(n): # O(log n)

    return fibonacci_n(n - 1) + fibonacci_n(n + 1)

def lucas(n):

    return (lucas_n(i) for i in xrange(n))

def catalan(n):

    return (1/(n + 1)) * scipy.misc.comb(2*n, n)

def bell(n):
    # Bell triange / Aitken's array / Peirce triangle
    '''
    :param n: natural number
    :return: list of n Bell Numbers
    '''
    triangle = np.zeros((n, n))
    triangle[0, 0] = 1
    for i in range(1, n):
        for j in range(i+1):
            if j == 0:
                triangle[i, j] = triangle[i - 1, i - 1]
            else:
                triangle[i, j] = triangle[i - 1, j - 1] + triangle[i, j - 1]

    return [triangle[k, k] for k in range(n)]

def bell_approx(n): # Dobinski's formula - nth moment of poisson distribution with E[x]=1
    '''
    :param n: natural number
    :return: nth Bell number
    '''

    mpmath.mp.dps = 50
    return int(round((1/math.e) * mpmath.nsum(lambda x: math.pow(x,n)/math.factorial(x), [0, mpmath.inf])))

def proth(a): # O(n^2)
    '''
    :param a: max block of proth numbers to return
    :return: 1-ath block of proth numbers
    '''

    pp = []
    pp.append(3)
    pp.append(5)
    count = 2
    increment = 3
    for i in range(1, a):
        for j in range(increment):
            pp.append(pow(2, i+1) + pp[count - 1])
            count += 1
        increment *= 2

    return pp



def sylvester(n): # O(n)
    '''
    :param n: number of sylvester numbers to return
    :return: n sylvester numbers
    '''
    num = []
    num.append(2)
    for i in xrange(1, n):
        num.append((num[i - 1] - 1) * num[i - 1] + 1)
    '''
    num = [(num[i - 1] - 1) * (num[i - 1] + 1) for i in xrange(1, n)]
    num.insert(0, 2)
    '''

    return num


def iscarmichael(n): # Carmichael numbers
    # knodel(1, m)
    if len(prime_factorization(n).keys()) >= 2: # Definition from Korselt
        p = divisors(n)
        div = [1 for i in p if (n-1) % (i-1) == 0]
        if np.prod(div) == len(p):
            return True
    return False



def is_knodel(n, m=10000): # Knodel numbers (subset of Carmichael numbers)
    # n: positive integer
    count1 = count2 = 0
    for i in xrange(1, m):
        if iscoprime(i, m):
            count1 += 1
            if pow(i, m - n, m) == 1:
                count2 += 1
    if count1 == count2:
        return True
    return False


def harmonic(n):
    # n > 0
    total = 0
    for i in xrange(1, n):
        total += 1/i
    return total

def partial_avg(s):
    # Given s != []
    items = []
    total = 0
    for i in xrange(len(s)):
        total += s[i]
        items.append(total / i)
    return items

def agm(a, g): # Arithmetic-geometric mean - real numbers
    # agm(1, sqrt(2)) -> Gauss's constant named after Carl Friedrich Gauss
    '''
    :param a: real number
    :param g: real number
    :return: arithmetic-geometric mean
    '''
    a1 = (a + g) / 2
    g1 = math.sqrt(a * g)
    while abs(a1 - g1) >= 0.0000000000000000001:
        an = (a1 + g1) / 2
        gn = math.sqrt(a1 * g1)
        a1 = an; g1 = gn

    return a1

def geom_harmonic(h, g): # Geometric-harmonic mean
    '''
    :param h: real number
    :param g: real number
    :return: geometric-harmonic mean
    '''
    h1 = 2/((1/h) + (1/g))
    g1 = math.sqrt(h * g)
    while abs(h1 - g1) >= 0.0000000000000000001:
        hn = 2/((1/h1) + (1/g1))
        gn = math.sqrt(h1 * g1)
        h1 = hn; g1 = gn

    return h1

def riemann_zeta(s):

    return mpmath.nsum(lambda x: 1./pow(x, s), [1, mpmath.inf])

def apery_const():

    return riemann_zeta(3)

from scipy.integrate import quad
def _integrand_catalan(x):

    return 0.5 * np.log(1/np.cos(x) + np.sin(x)/np.cos(x))

def catalan_const():

    return quad(_integrand_catalan, 0, math.pi/2)[0]

def dirichlet():
    pass

def pi(): # Gregory-Leibniz Series

    return 4 * mpmath.nsum(lambda n: pow(-1, n)/(2*n + 1), [0, mpmath.inf])

def pi_fast(): # Chudnovsky

    p = 12 * mpmath.nsum(lambda k: pow(-1, k)*math.factorial(6*k)*(13591409+545140134*k)/
                        (math.factorial(k*k)*pow(math.factorial(k),3)*pow(640320, 3*k+1.5)), [0, mpmath.inf])
    return 1/p

if __name__ == "__main__":
    print sylvester(5)
    print fibonacci_n(1000000)
    #print lucas_n(10000)
    #print iscarmichael(561)
    print riemann_zeta(2)
    print apery_const()
    print catalan_const()
    print "PI slow: " + str(pi())
    print "PI: " + str(pi_fast())