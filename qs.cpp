#include <gmpxx.h>
#include <vector>
#include <cmath>
#include <stdio.h>

// The optimal smoothness bound is exp((0.5 + o(1)) * sqrt(log(n)*log(log(n)))).
#define SMOOTH_BOUND 50
#define TRIAL_BOUND 100000
#define SIEVE_CHUNK 60


typedef std::vector<int> int_vector;
typedef std::vector<mpz_class> mpz_vector;


// printf due to some version incompatibility
void print_mpz_class(mpz_class x)
{
    gmp_printf("%Zd \n", x.get_mpz_t());
}

void print_mpz_vector(mpz_vector x)
{
    for(unsigned int i=0; i<x.size(); i++)
    {
        gmp_printf("%Zd, ", x[i].get_mpz_t());
    }
    gmp_printf("\n");
}

void print_int_vector(int_vector x)
{
    for(unsigned int i=0; i<x.size(); i++)
    {
        gmp_printf("%d, ", x[i]);
    }
    gmp_printf("\n");
}



// Sloppy coding
// Return a list of primes
int_vector eratosthenes(int bound)
{
    int_vector primes;

    std::vector<bool> A (bound, 1);

    A[0] = 0; // 0 and 1 aren't prime
    A[1] = 0;

    for(int i=2; i<sqrt(bound); i++)
    {
        if (A[i])
        {
            for(int j = i*i; j<=bound; j+=i)
            {
                A[j] = 0;
            }
        }
    }

    for(int i=0; i<bound; i++)
    {
        if (A[i])
        {
            primes.push_back(i);
        }
    }
    return primes;
}

// Return a vector of a number's factors (ex. [0, 1, 2, 0]) and a boolean of
// whether it's smooth or not
typedef std::pair<int_vector, bool> vb_pair;

vb_pair factor_smooth(mpz_class n, mpz_vector factor_base)
{
    // Each item in factors corresponds to number in factor base
    int_vector factors(factor_base.size(), 0);
    bool is_smooth = false;

    for(unsigned int i=0; i<factor_base.size(); i++)
    {
        mpz_class factor = factor_base[i];
        // print_mpz_class(factor);

        while (n % factor == 0)
        {
            n /= factor;
            factors[i] = (factors[i]+1) % 2; // mod 2 matrices
        }
    }
    if (n==1)
    {
        is_smooth = true;
    }
    vb_pair return_pair(factors, is_smooth);
    return return_pair;
}



int main()
{

    const mpz_class n = 90283;

    int_vector primes = eratosthenes(TRIAL_BOUND);

    mpz_vector factor_base;

    // Create factor base
    mpz_class two = 2;
    factor_base.push_back(two);
    for(unsigned int i=0; i<primes.size(); i++)
    {
        int p = primes[i];
        if (p > SMOOTH_BOUND) // Up to smooth limit
        {
            break;
        }
        mpz_class p_mpz = p;
        if (mpz_legendre(n.get_mpz_t(), p_mpz.get_mpz_t()) == 1)
        {
            //print_mpz_class(p);
            factor_base.push_back(p);
        }
    }


    // Find smooth numbers
    mpz_vector smooth_numbers;
    std::vector <int_vector> smooth_number_factors; // Corresponds to smooth numbers

    mpz_class j = 1; // x = sqrt(n) + j
    mpz_class sqrt_n = sqrt(n);

    //while (smooth_numbers.size() < factor_base + 1)
    //{
        mpz_vector current_chunk(SIEVE_CHUNK);
        for(int i=0; i<SIEVE_CHUNK; i++)
        {
            mpz_class current;
            mpz_class offset = sqrt_n+j+i; // Current addition to x
            // print_mpz_class(offset);
            // current = (j+i)^2 mod n
            mpz_powm_ui(current.get_mpz_t(), offset.get_mpz_t(), 2, n.get_mpz_t());
            //print_mpz_class(current);

            current_chunk[i] = current;

        }
        // To do: add Shanks-Tonelli
        j += SIEVE_CHUNK;
    //}

    //printf("%d", factor_base.size());

    //mpz_class test = 153;
    //vb_pair z = factor_smooth(test, factor_base);
    //gmp_printf("%d\n", z.second);


    for(unsigned int i=0; i<current_chunk.size(); i++)
    {
        vb_pair factored = factor_smooth(current_chunk[i], factor_base);
        if (factored.second) // Is smooth
        {
            smooth_numbers.push_back(current_chunk[i]);
            smooth_number_factors.push_back(factored.first);
        }
    }

    print_mpz_vector(smooth_numbers);
    for(int i=0; i<smooth_number_factors.size(); i++)
    {
        print_int_vector(smooth_number_factors[i]);
    }














    return 0;
}


