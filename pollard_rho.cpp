#include <iostream>
#include <gmpxx.h>

#define MIN(a, b) ((a < b) ? a : b)

// Euclid's algorithm
// Could also use GMP's implementation or binary GCD
inline mpz_class euclid_gcd(mpz_class a, mpz_class b)
{
    mpz_class c;
    while (a)
    {
        c = a; a = b % a; b = c;
    }
    return b;
}

inline mpz_class gmp_gcd(mpz_class a, mpz_class b)
{
    mpz_class c;
    mpz_gcd(c.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return c;
}

inline mpz_class func(mpz_class x, mpz_class c, mpz_class N)
{
    return ((x*x)%N + c) % N; // Using function f(x) = x*x + c
}

gmp_randclass r(gmp_randinit_default); // Initialize RNG



// Pollard's rho, based on Come On Code On's python code
mpz_class pollard_rho(mpz_class N)
{
    if (N % 2 == 0)
    {
        mpz_class two = 2;
        return two;
    }
    r.seed(123);

    mpz_class x = r.get_z_range(N-1) + 1; // 0 <= x <= N-1
    mpz_class y = x;
    mpz_class c = r.get_z_range(N-1) + 1;
    mpz_class g = 1;

    while (g==1)
    {
        x = func(x, c, N);
        y = func(y, c, N);
        y = func(y, c, N);
        g = gmp_gcd(abs(x-y), N);
    }
    return g;
}



// Brent's modification
mpz_class brent(mpz_class N)
{
    if (N % 2 == 0)
    {
        mpz_class two = 2;
        return two;
    }
    r.seed(123);

    mpz_class y = r.get_z_range(N-1) + 1;
    mpz_class c = r.get_z_range(N-1) + 1;
    mpz_class m = r.get_z_range(N-1) + 1;
    mpz_class g = 1, z = 1, q = 1, ys, x;

    while (g==1)
    {
        x = y;
        for(mpz_class i=0; i<z; i++)
        {
            y = func(y, c, N);
        }
        mpz_class k = 0;
        while (k<z && g==1)
        {
            ys = y;
            for(mpz_class i=0; i<MIN(m, z-k); i++)
            {
                y = func(y, c, N);
                q = q*(abs(x-y)) % N;
            }
            g = gmp_gcd(q, N);
            k += m;
        }
        z *= 2;
    }
    if (g==N)
    {
        while (true)
        {
            ys = func(ys, c, N);
            g = gmp_gcd(abs(x-ys), N);
            if (g>1)
                break;
        }
    }
    return g;
}

int main()
{
    const mpz_class N("10000000000009800000000002077");
    mpz_class p = pollard_rho(N);
    mpz_class q = N / p;
    std::cout << p << '\n' << q << '\n';

    return 0;
}
