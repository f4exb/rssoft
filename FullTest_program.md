# Full test program #

## Introduction ##

This test program does the following:

  1. Creates a random message of the specified length (-k option) corresponding to the k value in the RS(n,k) code
  1. Encodes the message into a codeword following RS(n,k) definition
  1. Optionnally punctures out a given number of erasures randomly (-e option)
  1. Applies Additive White Gaussian Noise (AWGN) to the remaining symbols and creates an analog "power" matrix with the n symbol positions as columns (x coordinate) and n+1 symbol values as rows (y coordinate). An erased symbol is represented as a column with all zeros. n+1=2^m is the size of the Galois Field where m is the log2 of the size (-m option)
  1. Normalizes the "power" matrix to obtain a reliability matrix
  1. Makes a multiplicity matrix according to given global multiplicity (-M option). The so called "global multiplicity" is the sum of all terms in the multiplicity mstrix
  1. Runs the Guruswami-Sudan with Koetter-Vardy variations interpolation algorithm using this multiplicity matrix to produce an interpolation bivariate polynomial
  1. Runs the Roth-Ruckenstein factorization algorithm to produce possible candidate encoding polynomials
  1. Estimate the codeword and calculate its a posteriori probability at all evaluation points excluding erasure locations. Order the results by decreasing probability
  1. Compare original message with corresponding encoding polynomial
  1. If it is found we are done
  1. Else return to multiplicity matrix calculation (step 6) taking the matrix cost as the new multiplicity matrix global multiplicity until the maximum number of retries (-i option) is reached
  1. Optionnally (--print-stats option) print the run data on a single line starting with "_RES" for scripted automatic runs._

## List of parameters ##

| option | descrption | default |
|:-------|:-----------|:--------|
| --print-seed | Outputs the value of the random seed being used | false |
| --print-stats | Outputs the run data on a single line starting with "#RES" for scripted automatic runs | false |
| --sagemath | Outputs the input parameters to the Sage decoding script | false |
| -n | Signal to noise ratio applied to the codeword samples | 0 |
| -m | Log2 of Galois Field size i.e. GF(q) with q = 2^m | 3 |
| -k | Number of message symbols as in RS(n,k) | 5 |
| -M | Global multiplicity i.e. sum of multiplicity of all points in matrix (x,y) space | 12 |
| -v | Verbosity level from 0 (minimal) to 2 (very verbose) | 0 |
| -s | If present seed for random engine. This ensures the same random conditions can be reprduced | (none) |
| -i | Maximum number of retries | 0 |
| -e | Number of erasures | 0 |

## Statistics print format ##

This is a single line starting with "@RES " with following parameters in a comma separated list in this order:
  1. 1 if codeword was found else 0
  1. Number of trials
  1. Probability score order of the result codeword at the successful attempt
  1. Number of result codewords obtained at the successful attempt
  1. Codeword average probability score in dB/symbol
  1. Number of hard decision errors
  1. Number of erasures specified
  1. Number of false results
  1. SNR specified
  1. Maximum matrix global multiplicity
  1. Maximum multiplicity matrix cost

## Example ##

```
FullTest -m 3 -k 3 -n 4.7 -M 8 -i 5 -e 3 --print-seed --print-stats
Seed = 1731716427
use seed: 1731716427
Message : (k=3) [2, 3, 3]
Codeword: (n=7) [2, 3, 4, 3, 5, 5, 4]
Erasures: (n=7) [*, *, *, 3, 5, 5, 4]
Hard-dec: (n=7) [*, *, *, 3, 5, 5, 3]
 -> 1 errors, 3 erasures: uncorrectable with hard decision
Codeword score: -5.27202 dB/symbol (best = -2.38352, worst = -11.1738)

Multiplicity matrix cost is 10
Q(X,Y) = a^5 + a^1*X + a^5*Y + a^4*X^2 + a^1*X*Y + a^3*X^3 + a^4*Y^2 + a^2*X^2*Y + a^2*X^4 + a^5*X^3*Y + a^6*X^5
0 result(s)

Multiplicity matrix cost is 13
Q(X,Y) = a^2 + X + a^1*Y + a^3*X^2 + a^4*X*Y + a^2*X^3 + a^1*X^2*Y + a^4*X^4 + a^2*X^3*Y + a^4*X^5 + a^4*X^4*Y + X^6
1 result(s)
F0(X) = a + a^3*X + a^3*X^2
Codewords:
#0: (-5.27 dB/symbol) [2, 3, 4, 3, 5, 5, 4]
Messages:
#0: (-5.27 dB/symbol) [2, 3, 3]
#0 found at iteration #2 !!!

#RES: 4.7,-5.27202,1,3,1,1,0,2,0,10,13
```