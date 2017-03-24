/* 
Uses mean parameterization, i.e. tau = 1/rate.

*/

/** \brief Vectorize 5-argument functions.

    For five-argument functions (Type, Type, Type, Type, int),
    vectorize first four arguments.
*/
#define VECTORIZE5_tttti(FUN) \
  GVECTORIZE(FUN,V,T,T,T,I,N) \
  GVECTORIZE(FUN,T,V,T,T,I,N) \
  GVECTORIZE(FUN,T,T,V,T,I,N) \
  GVECTORIZE(FUN,T,T,T,V,I,N) \
  GVECTORIZE(FUN,V,V,T,T,I,N) \
  GVECTORIZE(FUN,T,V,V,T,I,N) \
  GVECTORIZE(FUN,T,T,V,V,I,N) \
  GVECTORIZE(FUN,V,V,V,T,I,N) \
  GVECTORIZE(FUN,T,V,V,V,I,N) \
  GVECTORIZE(FUN,V,V,V,V,I,N) \
  GVECTORIZE(FUN,V,T,T,V,I,N)

#define VECTORIZE3_n(FUN) \
template<class Type> \
vector<Type> FUN(int n, Type arg1, Type arg2, Type arg3) { \
    vector<Type> ans(n); \
    for(int i=0; i<n; i++) ans(i) = FUN(arg1, arg2, arg3); \
    return ans; \
}

namespace exgaussian {
    template<class Type>
    Type dexgaussian(Type x, Type mu, Type sigma,
                     Type tau, int give_log = 0)
    {
        Type temp;
        Type out;
        temp = x - mu - (pow(sigma, 2)/tau);
        out = -log(tau) - (temp + (pow(sigma, 2)/(2 * tau)))/tau + 
                   log(pnorm(temp/sigma));
        if (!give_log) 
            return CppAD::CondExpGe(x, Type(0), exp(out), Type(0));
        else 
            return CppAD::CondExpGe(x, Type(0), out, Type(-INFINITY));
        
    }
    VECTORIZE5_tttti(dexgaussian)

    template<class Type>
    Type pexgaussian(Type q, Type mu, Type sigma,
                     Type tau, int lower_tail = 1, int give_log = 0)
    {
        Type sig_sq = pow(sigma, 2);
        Type z = q - mu - sig_sq/tau;
        Type out = pnorm((q - mu)/sigma) - pnorm(z/sigma) *
              exp(pow((mu + sig_sq)/tau, 2) - pow(mu, 2) - 2 * q * sig_sq/tau) / (2 * sig_sq));
        if (!lower_tail)
            out = 1 - out;
        if (give_log)
            out = log(out);
        return(out);
    }

    template<class Type>
    Type rexgaussian(Type mu, Type sigma, Type tau) 
    {
        return (rnorm(mu, sigma) + rexp(tau));
    }
    VECTORIZE3_n(rexgaussian)
}