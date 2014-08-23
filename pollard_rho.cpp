#include <algorithm>
#include <gmpxx.h>

#define MIN(a, b) ((a < b) ? a : b)

// Euclid's algorithm
// Based on Come On Code On's python code
inline mpz_class gcd(mpz_class a, mpz_class b)
{
    mpz_class c;
    while (a)
    {
        c = a; a = b % a; b = c;
    }
    return b;
}
// Initialize random number generator
gmp_randclass r(gmp_randinit_default);


// Pollard's rho using f(x) = x*x + c
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
        x = ((x*x)%N + c) % N;
        y = ((y*y)%N + c) % N;
        y = ((y*y)%N + c) % N;
        g = gcd(abs(x-y), N);
    }
    return g;
}

// Brent's modification
mpz_class pollard_brent_rho(mpz_class N)
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
            y = ((y*y)%N + c) % N;
        }
        mpz_class k = 0;
        while (k<z && g==1)
        {
            ys = y;
            for(mpz_class i=0; i<MIN(m, z-k); i++)
            {
                y = ((y*y)%N + c) % N;
                q = q*(abs(x-y)) % N;
            }
            g = gcd(q, N);
            k += m;
        }
        z *= 2;
    }
    if (g==N)
    {
        gmp_printf("hi \n");
        while (true)
        {
            ys = ((ys*ys)%N + c) % N;
            g = gcd(abs(x-ys), N);
            if (g>1)
                break;

        }
    }
    return g;
}

int main()
{
    const mpz_class N("519920418481914614701");
    mpz_class p = pollard_brent_rho(N);
    mpz_class q = N / p;
    gmp_printf("%Zd\n%Zd", p.get_mpz_t(), q.get_mpz_t());

    return 0;
}
