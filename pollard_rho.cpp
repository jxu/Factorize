#include <iostream>
#include <gmpxx.h>
#include <ctime>

#define MIN(a, b) ((a < b) ? a : b)

void *__gxx_personality_v0;
void *_Unwind_Resume;

// Euclid's algorithm (much slower than GMP's built-in)
inline mpz_class euclid_gcd(mpz_class a, mpz_class b)
{
    mpz_class c;
    while (a)
    {
        c = a; a = b % a; b = c;
    }
    return b;
}


inline void func(mpz_class &x, const mpz_class &c, const mpz_class &N)
{
    x = (x*x + c) % N; // Using function f(x) = x*x + c
}

gmp_randclass r(gmp_randinit_default); // Initialize RNG


// Pollard's rho, based on Come On Code On's python code and activestate
mpz_class pollard_rho(const mpz_class &N)
{
    if (N % 2 == 0)
    {
        mpz_class two = 2;
        return two;
    }

    //mpz_class x = r.get_z_range(N-1) + 1; // 0 <= x <= N-1
    mpz_class x = 2;
    mpz_class y = x;
    mpz_class c = r.get_z_range(N-1) + 1;
    mpz_class g = 1;
    mpz_class a; // temp variable

    std::cout << "c:" << c << std::endl;

    while (g==1)
    {
        func(x, c, N);
        func(y, c, N);
        func(y, c, N);
        a = abs(x-y);
        mpz_gcd(g.get_mpz_t(), a.get_mpz_t(), N.get_mpz_t());
    }
    return g;
}


// Brent's modification
mpz_class brent(const mpz_class &N)
{
    if (N % 2 == 0)
    {
        mpz_class two = 2;
        return two;
    }

    mpz_class y = r.get_z_range(N-1) + 1;
    mpz_class c = r.get_z_range(N-1) + 1;
    mpz_class m = r.get_z_range(N-1) + 1;
    mpz_class g = 1, z = 1, q = 1;
    mpz_class ys, x, k, a, i;

    std::cout << "y:" << y << "\nc:" << c << "\nm:" << m << std::endl;

    while (g==1)
    {
        x = y;
        for(i=0; i<z; i++)
            func(y, c, N);

        k = 0;
        while (k<z && g==1)
        {
            ys = y;
            for(i=0; i<MIN(m, z-k); i++)
            {
                func(y, c, N);
                q = (q * abs(x-y)) % N;
            }
            mpz_gcd(g.get_mpz_t(), q.get_mpz_t(), N.get_mpz_t());
            //g = gmp_gcd(q, N);
            k += m;
        }
        z *= 2;
    }
    if (g==N)
    {
        do
        {
            func(ys, c, N);
            a = abs(x-ys);
            mpz_gcd(g.get_mpz_t(), a.get_mpz_t(), N.get_mpz_t());
        } while (g == 1);
    }
    return g;
}

int main()
{
    r.seed(long(std::time(0)));
    const mpz_class N("10000000000009800000000002077");
    mpz_class p = brent(N);
    mpz_class q = N / p;
    std::cout << p << '\n' << q << '\n';

    return 0;
}
