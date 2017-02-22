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
    nll -= exgaussian::dexgaussian(y, mu, sigma, nu, true).sum();

    SIMULATE {
        y = exgaussian::rexgaussian(y.size(), mu, sigma, nu);
        REPORT(y);
    }
    return nll;
}
