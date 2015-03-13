# A few observations #

The solution designed and implemented here generally works and can routinely achieve results beyond the theoretical hard decision limit.

However while trying out on various examples the following observations can be made:

## Interpolation solution ##

While the interpolation algorithm is always successful in finding a bivariate polynomial solution this solution does not necessarily yields roots in Y and therefore there are no solutions to the subsequent factorization algorithm.

Sometimes the bivariate polynomial has no terms in Y effectively making it an univariate polynomial in X which obviously cannot have any root in Y. This has been observed only on short codes (i.e. over GF(8) field).

More commonly the bivariate polynomial which has Y coefficients has no root in Y. This would happen for any GF size but only in difficult conditions.

## Multiplicity ##

In the following the term "global multiplicity" designates the sum of all elements in the multiplicity matrix.

When no solution to the factorization (i.e. no Y root) is found at a given global multiplicity it may happen that a solution exist and can be successful for a larger value of the global multiplicity. Empirically one would try again by taking the multiplicity matrix cost as the next global multiplicity.

The improvement in finding solution as the global multiplicity increases is in no way linear. When reaching a value that yields solutions sometimes a larger value for the global multiplicity will not give any solution.

## Other caveats ##

If the message is all zero symbols then the encoded codeword is also all zeros and there is no interpolation solution. So the message with all zeros cannot be used.