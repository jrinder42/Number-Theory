import numpy as np
from Primes import isprime, prime_factorization, divisors
import math
import cmath



def legendre_symbol(a, p):
    '''
    :param a: integer
    :param p: positive prime
    :return: integer -1, 0, or 1
    '''

    if isprime(p) == False or p <= 2:
        raise ValueError("Incorrect number passed")

    quad_residue = list(set([pow(n, 2) % p for n in range(0, int(math.floor(p/2.)+1))]))

    if a % p == 0:
        return 0
    elif a % p in quad_residue:
        return 1
    else:
        return -1

def jacobi_symbol(a, n): # Generalization of legendre symbol
    '''
    :param a: integer
    :param n: positive integer
    :return: integer -1, 0, or 1
    '''

    if n % 2 == 0:
        return ValueError

    primes = prime_factorization(n).keys()
    powers = prime_factorization(n).values()
    return np.prod([pow(legendre_symbol(a, primes[i]), powers[i]) for i in range(len(primes))])

def kronecker_symbol(a, n): # Generalization of jacobi symbol
    '''
    :param a - integer:
    :param n - nonzero integer:
    :return - integer -1, 0, or 1:
    '''

    if n == 0:
        return ValueError
    pre = []
    u = 1
    if n < 0:
        u = -1

    if n == 0:
        if a == 1 or a == -1:
            pre.append(1)
        else:
            pre.append(0)

    if u == 1 and n == 1:
        pre.append(1)

    if u == -1 and n == -1:
        if a < 0:
            pre.append(-1)
        else:
            pre.append(1)

    if (n > 1 or n < 1) and n != 0:
        factors = divisors(n)
        for p in factors:
            if p % 2 != 0:
                pre.append(legendre_symbol(a, p))
            else:
                if a % 2 == 0:
                    pre.append(0)
                elif a % 8 == 1 or a % 8 == 7:
                    pre.append(1)
                elif a % 8 == 3 or a % 8 == 5:
                    pre.append(-1)

    return np.prod(pre)

def cipolla(a, p):
        '''
        :param a - integer:
        :param p - positive odd prime:
        :return - list containing cipolla values:
        '''

        if legendre_symbol(a, p) == 1 and isprime(p):
            Fp = range(p)
            a = 0
            for i in Fp:
                print legendre_symbol(i**2 - a, p)
                if legendre_symbol(i**2 - a, p) == -1:
                    a = i
                    break
            modify = pow(a + cmath.sqrt(a**2 - a), (p+1)/2)
            alpha = int(round(modify.real % p))
            # beta = int(round(modify.imag / cmath.sqrt(a**2 - n).imag)) % p
            return [alpha, p - alpha]
        else:
            return ValueError

if __name__ == "__main__":
    pass