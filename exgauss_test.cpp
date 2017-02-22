#include <TMB.hpp>
#include "exgaussian.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
    DATA_VECTOR(y);
    PARAMETER(mu);
    PARAMETER(sigma);
    PARAMETER(nu);

    Type nll = 0;
    for(int total = 0; total < y.size(); total++) {
        nll -= exgaussian::dexgaussian(y[total], mu, sigma, nu, true);
    }

    SIMULATE {
        for(int total = 0; total < y.size(); total++) {
            y[total] = exgaussian::rexgaussian(mu, sigma, nu);
        }
        REPORT(y);
    }
    return nll;
}
