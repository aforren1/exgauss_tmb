/* 
Uses mean parameterization, i.e. nu = 1/rate.

*/

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

    template<class Type>
    Type rexgaussian(Type mu, Type sigma, Type nu) 
    {
        return (rnorm(mu, sigma) + rexp(1/nu));
    }
}