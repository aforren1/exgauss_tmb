/* 
Uses mean parameterization, i.e. nu = 1/rate.

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
                     Type nu, int give_log = 0)
    {
        Type temp;
        Type out;
        temp = x - mu - (pow(sigma, 2)/nu);
        out = -log(nu) - (temp + (pow(sigma, 2)/(2 * nu)))/nu + 
                   log(pnorm(temp/sigma));
        if (give_log) return out; else return exp(out);
    }
    VECTORIZE5_tttti(dexgaussian)

    template<class Type>
    Type rexgaussian(Type mu, Type sigma, Type nu) 
    {
        return (rnorm(mu, sigma) + rexp(nu));
    }
    VECTORIZE3_n(rexgaussian)
}