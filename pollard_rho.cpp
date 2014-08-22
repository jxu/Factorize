#include <iostream>
#include <gmpxx.h>

// Euclid's algorithm
inline mpz_class gcd(mpz_class a, mpz_class b)
{
    mpz_class c;
    while (a)
    {
        c = a; a = b % a; b = c;
    }
    return b;
}



mpz_class pollard_rho(mpz_class N)
{
    if (N % 2 == 0)
    {
        mpz_class two = 2;
        return two;
    }
    gmp_randclass r(gmp_randinit_default);
    r.seed(0);
    mpz_class x = r.get_z_range(N-1) + 1; // 1 <= x <= N-1
    mpz_class y = x;
    mpz_class c = r.get_z_range(N-1) + 1;
    mpz_class g = 1;

    while (g==1)
    {
        // f(x) = x*x + c
        x = ((x*x)%N + c)%N;
        y = ((y*y)%N + c)%N;
        y = ((y*y)%N + c)%N;
        g = gcd(abs(x-y), N);
    }
    return g;
}

int main()
{
    const mpz_class N("519920418481914614701");
    mpz_class p = pollard_rho(N);
    mpz_class q = N / p;
    gmp_printf("%Zd\n%Zd", p.get_mpz_t(), q.get_mpz_t());

    return 0;
}
