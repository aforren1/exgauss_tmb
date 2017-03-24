#include <TMB.hpp>
#include "exgaussian.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
    DATA_VECTOR(y);
    PARAMETER(mu);
    PARAMETER(sigma);
    PARAMETER(tau);

    Type nll = 0;
    nll -= exgaussian::dexgaussian(y, mu, sigma, tau, true).sum();

    SIMULATE {
        y = exgaussian::rexgaussian(y.size(), mu, sigma, tau);
        REPORT(y);
    }
    return nll;
}
