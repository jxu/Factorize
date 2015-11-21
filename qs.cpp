#include <gmpxx.h>
#include <vector>
#include <cmath>
#include <iostream>

// Constants
// The optimal smoothness bound is exp((0.5 + o(1)) * sqrt(log(n)*log(log(n)))).
const int SMOOTH_BOUND = 500;
const int TRIAL_BOUND = 400;
const int SIEVE_CHUNK = 60;

const bool DEBUG = true;



typedef std::vector<int> int_vector;
typedef std::vector<int_vector> matrix;
typedef std::vector<mpz_class> mpz_vector;

template <typename T> // Takes int_vector or mpz_vector
void print_vector(const T &x)
{
    for(size_t i=0; i<x.size(); i++)
        std::cout << x[i] << ", ";

    std::cout << '\n';
}


// Sloppy coding
// Return a list of primes
int_vector eratosthenes(int bound)
{
    int_vector primes;

    std::vector<bool> A (bound, 1);
    A[0] = 0; A[1] = 0; // 0 and 1 aren't prime

    for(int i=2; i<sqrt(bound); i++)
    {
        if (A[i])
        {
            for(int j = i*i; j<=bound; j+=i)
                A[j] = 0;
        }
    }

    for(int i=0; i<bound; i++)
    {
        if (A[i])
            primes.push_back(i);
    }
    return primes;
}

// Return a vector of a number's factors (ex. [0, 1, 2, 0]) and a boolean of
// whether it's smooth or not
typedef std::pair<int_vector, bool> vb_pair;
vb_pair factor_smooth(mpz_class n, const mpz_vector &factor_base)
{
    // Each item in factors corresponds to number in factor base
    int_vector factors(factor_base.size(), 0);

    for(size_t i=0; i<factor_base.size(); i++)
    {
        mpz_class factor = factor_base[i];
        while (n % factor == 0)
        {
            n /= factor;
            factors[i] ^= 1; // + 1 (mod 2) matrices
        }
    }
    bool is_smooth = (n==1);
    vb_pair return_pair(factors, is_smooth);
    return return_pair;
}

// ------------------------------------------------

int main()
{
    // Test numbers: 502560280658509, 90283
    const mpz_class n("502560280658509");

    int_vector primes = eratosthenes(TRIAL_BOUND);
    mpz_vector factor_base;

    // Create factor base
    mpz_class two = 2;
    factor_base.push_back(two);
    for(size_t i=0; i<primes.size(); i++)
    {
        int p = primes[i];
        if (p > SMOOTH_BOUND) // Up to smooth limit
            break;
        mpz_class p_mpz = p;
        // Use only primes that match (n|p) = 1
        if (mpz_legendre(n.get_mpz_t(), p_mpz.get_mpz_t()) == 1)
        {
            factor_base.push_back(p);
        }
    }

    if (DEBUG)
    {
       std::cout << "Factor base: ";
       print_vector(factor_base);
    }


    // Find smooth numbers with x = sqrt(n) + j
    mpz_class j = 1;
    mpz_class sqrt_n = sqrt(n);

    // Allocate factor_base+1 size
    mpz_vector smooth_numbers(factor_base.size()+1, 0);
    mpz_vector smooth_x(factor_base.size()+1, 0);

    // Corresponds to smooth numbers
    matrix smooth_factors(factor_base.size()+1, {0});
    int smooth_count = 0;
    bool not_done = true;

    while (not_done) // pi(B) + 1
    {
        mpz_vector current_chunk(SIEVE_CHUNK);
        mpz_vector current_x(SIEVE_CHUNK);
        for(int i=0; i<SIEVE_CHUNK; i++)
        {
            mpz_class current;
            mpz_class x = sqrt_n + j + i; // Current addition to x
            // current = (j+i)^2 mod n
            mpz_powm_ui(current.get_mpz_t(), x.get_mpz_t(), 2, n.get_mpz_t());

            current_chunk[i] = current;
            current_x[i] = x;
        }
        // To do: add Shanks-Tonelli
        j += SIEVE_CHUNK;

        // Actual factoring
        for(size_t i=0; i<current_chunk.size(); i++)
        {
            vb_pair factored = factor_smooth(current_chunk[i], factor_base);
            if (factored.second) // Is smooth
            {
                if (smooth_count > factor_base.size())
                {
                    not_done = false;
                    break;
                }
                smooth_x[smooth_count] = current_x[i];
                smooth_numbers[smooth_count] = current_chunk[i];
                smooth_factors[smooth_count] = factored.first;
                smooth_count++;
            }
        }
    }


    if (DEBUG)
    {
        std::cout << "Smooth x: ";
        print_vector(smooth_x);
        std::cout << "Smooth numbers: ";
        print_vector(smooth_numbers);

        std::cout << "Smooth factors:\n";
        for(size_t i=0; i<smooth_factors.size(); i++)
            print_vector(smooth_factors[i]);

        std::cout << '\n';
    }


    // Gaussian Elimination -----------------------------------
    // Transpose the matrix
    int Ai = smooth_factors[0].size(); // row
    int Aj = smooth_factors.size(); // column
    matrix A(Ai, int_vector(Aj, 0));

    for(int i=0; i<Ai; i++)
        {
            for(int j=0; j<Aj; j++)
            {
                A[i][j] = smooth_factors[j][i];
            }
        }

    if (DEBUG)
    {
        std::cout << "Transposed matrix A:\n";
        for(size_t i=0; i<A.size(); i++)
            print_vector(A[i]);

            std::cout << '\n';
    }

    for(int k=0; k<Ai; k++)
    {
        // Swap with pivot if current diagonal is 0
        if (A[k][k] == 0)
        {
            for(int l=k; l<Ai; l++)
            {
                if (A[l][k]==1)
                {
                    A[l].swap(A[k]);
                    break;
                }
            }
        }
        // For rows below pivot
        for(int i=k+1; i<Ai; i++)
        {
            // If row can be subtracted, subtract every element (using xor)
            if (A[i][k])
            {
                for(int j=0; j<Aj; j++)
                    A[i][j] ^= A[k][j];
                //for(size_t i=0; i<A.size(); i++)
                //   print_vector(A[i]);
                //std::cout << '\n';
            }
        }
    }

    // Find line between free and pivot variables
    int f;
    for(f=0; f<Aj; f++)
    {
        if (A[f][f] != 1)
            break;
    }

    // Back substitution on upper triangular matrix
    for(int k=f-1; k>=0; k--)
    {
        for(int i=k-1; i>=0; i--)
        {
            if (A[i][k])
            {
                for(int j=0; j<Aj; j++)
                    A[i][j] ^= A[k][j];
            }
        }
    }

    if (DEBUG)
    {
        std::cout << "Fully reduced matrix:\n";
        for(size_t i=0; i<A.size(); i++)
            print_vector(A[i]);

        std::cout << '\n';
    }


    int_vector null_space(Aj, 0);
    // Subject to change
    // First free variable is 1, rest are 0
    null_space[f] = 1;

    for(int i=0; i<f; i++)
            null_space[i] = A[i][f];


    mpz_class x_square = 1;
    mpz_class y = 1;
    for(size_t i=0; i<null_space.size(); i++)
        if (null_space[i])
        {
            x_square *= smooth_numbers[i];
            y *= smooth_x[i];
        }

    if (DEBUG)
    {
        std::cout << "Null space: ";
        print_vector(null_space);
        std::cout << "Square: " << x_square << std::endl;

    }
    mpz_class x, rem;
    mpz_sqrtrem(x.get_mpz_t(), rem.get_mpz_t(), x_square.get_mpz_t());


    if (DEBUG)
    {
        if (rem==0)
            std::cout << "Remainder 0\n";

        std::cout << "x: " << x << '\n' << "y: " << y << "\n\n";
    }


    mpz_class factor_1;
    mpz_class dif = y - x;
    mpz_gcd(factor_1.get_mpz_t(), n.get_mpz_t(), dif.get_mpz_t());
    if (factor_1 == 1 || factor_1 == n)
        std::cout << "Factoring failure: try again with different parameters\n";

    mpz_class factor_2 = n / factor_1;

    std::cout << "Factor 1: " << factor_1 << '\n';
    std::cout << "Factor 2: " << factor_2 << '\n';

    return 0;
}
