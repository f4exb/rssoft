# Full test program for libccsoft testing #

## Introduction ##

This test program does the following:

  1. Generates the required number of input symbols (-r option) or uses the specified list of input symbols (-i option). The longest constraint length (longest register length minus one) of zeroes is appended
  1. Encodes the message into a (n,k,m) convolutional code
  1. Applies Additive White Gaussian Noise (AWGN) to the obtained symbols and creates an analog "power" matrix with the output symbol positions as columns (x coordinate) and output symbol values (n bit symbols) as rows (y coordinate).
  1. Normalizes the "power" matrix to obtain a reliability matrix
  1. Runs the convolutional decoding algorithm using the symbol reliability information from the reliability matrix
  1. Compare original message with input symbols corresponding to the path elected by the algorithm
  1. Along with documented information, print the run data on a single line starting with "`_`RES" for scripted automatic runs.

With the Fano algorithm implementation there is an optional code tree cache mechanism that works basically as follows:

  * Nodes (and edges) are never deleted when they are abandoned (move backwards)
  * Before creating new nodes the node limit in tree cache is tested
    * If it is below limit the child nodes and edges are allocated in memory
    * If it is above the limit a purge mechanism is engaged before allocating new nodes and edges. This purge mechanism keeps only the current path and all path nodes immediate successors. This results in the same tree as if the abandoned nodes deletion was in place (no cache).

## Variants ##

  * **FullTest** uses std::vectors for registers and forward node+edge combos pointers.

  * **FullTest\_FA** uses std::arrays for registers and forward node+edge combos pointers. It uses the classes of CCSoft with _FA suffix._

## List of parameters ##

| option | descrption | default |
|:-------|:-----------|:--------|
| --print-seed | Outputs the value of the random seed being used | false |
| -n | Signal to noise ratio applied to the codeword samples | 0 |
| -v | Verbosity level from 0 (minimal) to 2 (very verbose) | 0 |
| -d | Specify the file name of a Graphviz dot file representing the explored code tree | (none) |
| -k | Comma separated list of constraint lengths (actually register lengths) for each of the k input lines | (none) |
| -g | Colon separated lists (comma separated) of generator polynomials for each of the n output lines. Each list has as many elements as there are input lines (k).| (none) |
| -i | Comma separated list of input symbols | (none) |
| -r | Generate this number of random symbols  | 0 |
| -s | Specify random generator seed | (not specified) |
| -N | Specify number of created nodes limit | (not specified) |
| -M | Specify a lower metric limit | (not specified) |
| -a | Specify the decoding algorithm and its parameters (see next) | stack:-1 |

## Decoding algorithm specification ##
### Stack algorithm ###
Format is `-a stack:<metric bias>`
The edge metric is defined as the base 2 logarithm of the edge's output symbol reliability obtained from the reliability matrix plus the metric bias
### Fano algorithm ###
Format is `-a fano:<metric bias>,<initial metric threshold>,<metric delta>,<tree cache size>,<initial metric threshold delta>`
  * metric bias: Defined in the same way as for the stack algorithm
  * initial metric threshold: The algorithm starts with this metric threshold. The optimal value has to be searched experimentally between a too low value where the algorithm cannot reach the end of the tree and a too high limit where the loop condition is detected. A loop condition is detected when the algorithm reaches back the root node in the same conditions as the starting conditions.
  * metric delta: This is the amount by which the threshold is tightened or loosened. A commonly working value is 1.
  * tree cache size: if given and positive this is the maximum number of nodes that can be cached (0 means no memory cache)
  * initial metric threshold delta: if specified when a loop condition is detected instead of aborting it restarts the whole process from the root node with an initial metric threshold added with this (negative) amount (hence decremented). If amount is positive or zero then the feature is not activated.

## Example ##
### Layland and Lushbaugh rate 1/2 with 32 bit register ###
This is the code used in the [JT4 protocol](http://physics.princeton.edu/pulsar/K1JT/JT2_JT4.TXT) for weak signal communications. We will use the Fano algorithm because of the large number of symbols.
#### Without tree cache ####
  * code rate is 1/2 with 1 bit input symbols and 2 bits output symbols
  * constraint length is 31 (32 bit register)
  * generator polynomials are 0xf2d05351 and 0xe4613c47 which corresponds to decimal values of 4073739089 and 3831577671 respectively. This defines the Layland and Lushbaugh code.
  * we use a metric bias of -1. For 2 bit output symbols (4 values) the equiprobability is 1/4 for which log 2 is -2. A value of -1 corresponds to the log 2 of twice this equiprobability value.
  * We use a starting metric threshold value of -20 which seems to give good results
  * We use a threshold delta of 1
  * we generate 72 input symbols as in the JT type source coded message to which 31 zero symbols are added
  * AWGN noise is added directly to the soft symbol values for a SNR of 3.4 dB
```
FullTest -k 32 -g 4073739089,3831577671 -r 72 -n 3.4 -N 10000000 -a fano:-1,-20,1
k=1, n=2, m=32
0 (32) : f2d05351 e4613c47 
1 1 1 1 1 1 0 0 0 1 1 0 0 0 0 1 0 0 0 0 1 0 0 1 0 1 1 1 1 1 0 0 0 1 1 0 1 1 1 0 1 0 1 1 1 0 1 1 1 1 1 1 0 1 1 0 1 1 0 0 0 0 1 0 1 1 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
3 1 3 3 2 2 2 0 3 1 0 3 1 0 0 0 1 1 1 1 0 1 3 0 2 1 1 0 3 1 2 0 0 3 3 0 2 0 1 1 1 3 1 3 3 3 3 2 3 0 2 0 1 2 3 3 3 0 3 1 3 3 3 1 2 3 1 0 2 2 2 0 1 2 0 1 0 2 1 0 1 1 1 1 1 2 1 0 2 2 3 0 1 0 1 1 1 1 3 1 3 3 3 
[1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] Success!
score = -4.9434 cur.threshold = -20 nodes = 1249606 eff.nodes = 176 moves = 1447998 max depth = 102
_RES 1,-4.9434,-20,1249606,176,1447998,102
```

#### With tree cache ####
  * we use a limit of 200M node creations
  * the tree cache has a limit of 4M nodes
  * it is successful at last after 66M node creations and 75M moves!
  * this is a debug version displaying the purge of the code tree cache each time it reaches the limit of 4M nodes.
  * htop running in parallel shows a maximum resident memory of 743 MB. This gives a cache size of about **185 MB per million of nodes**.
```
FullTest -k 32 -g 4073739089,3831577671 -r 72 -n 3.4 -N 200000000 -a fano:-1,-20,1,4000000
k=1, n=2, m=32
0 (32) : f2d05351 e4613c47 
1 0 0 1 1 0 0 0 1 0 0 1 1 1 0 0 1 1 0 1 1 1 1 0 1 0 1 1 1 1 0 0 1 1 0 0 0 1 1 0 0 1 1 1 1 0 1 1 1 1 0 1 0 0 1 0 0 1 1 0 0 0 0 1 0 0 1 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
3 2 2 3 0 0 1 1 3 0 3 0 3 2 2 2 0 1 0 3 0 0 3 1 3 0 3 1 3 0 3 2 2 3 2 3 3 1 2 1 3 2 3 1 3 1 2 1 0 3 3 1 3 0 2 3 0 3 2 1 0 1 1 2 3 3 0 1 3 3 1 1 2 2 1 0 1 3 1 0 2 1 1 2 0 3 0 2 3 3 0 0 3 3 1 0 0 1 1 3 0 3 0 
purged tree cache, nb of remaining nodes = 135
purged tree cache, nb of remaining nodes = 121
purged tree cache, nb of remaining nodes = 146
purged tree cache, nb of remaining nodes = 137
purged tree cache, nb of remaining nodes = 143
purged tree cache, nb of remaining nodes = 137
purged tree cache, nb of remaining nodes = 125
purged tree cache, nb of remaining nodes = 137
purged tree cache, nb of remaining nodes = 143
purged tree cache, nb of remaining nodes = 143
purged tree cache, nb of remaining nodes = 117
purged tree cache, nb of remaining nodes = 121
purged tree cache, nb of remaining nodes = 139
purged tree cache, nb of remaining nodes = 123
purged tree cache, nb of remaining nodes = 99
purged tree cache, nb of remaining nodes = 129
[1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] Success!
score = 9.55823 cur.threshold = -20 nodes = 66175332 eff.nodes = 2177419 moves = 75059495 max depth = 102
_RES 1,9.55823,-20,66175332,2177419,75059495,102
```

#### With unloop feature ####

When a loop condition is detected, that is the process returns to root node with the same threshold limit as at the beginning, instead of aborting it restarts with a lower threshold. The lowering amount is given as a negative number, here -4.

The metric limit (-M option) is given to limit the number of such retries. Below a certain limit there is no much hope to find a solution in a decent time. A node number threshold limit (-N option) can also be given to limit the process time.

This "unloop" or retry feature can be a good compromise between finding a solution quickly in favorable conditions without giving up too easily on difficult cases. It has been proven to be the best compromise between processing time and sensitivity.

```
FullTest -k 32 -g 4073739089,3831577671 -r 72 -n 3.6 -N 100000000 -a fano:-1,-8,1,10000000,-4 -M -24
k=1, n=2, m=32
0 (32) : f2d05351 e4613c47 
0 1 0 1 1 0 1 1 1 1 0 1 0 1 1 1 0 0 1 0 0 1 0 1 0 1 0 1 0 1 0 0 1 0 0 0 1 1 0 0 0 1 1 1 1 0 1 1 0 0 1 0 1 1 0 0 1 1 1 1 0 1 1 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 3 2 1 1 1 1 3 2 1 3 3 2 3 1 1 3 1 1 2 0 2 2 1 2 1 2 0 2 0 3 3 0 1 1 1 3 0 0 0 0 3 1 3 3 3 2 0 0 3 3 3 1 0 3 1 3 1 0 3 0 0 2 0 2 3 0 1 3 1 0 1 2 0 2 1 1 2 3 3 3 0 0 1 3 2 1 0 3 2 3 0 2 3 1 1 2 0 0 3 0 0 0 
Loop condition detected, restart with init threshold = -12
Loop condition detected, restart with init threshold = -16
[0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] Success!
score = -10.7666 cur.threshold = -16 nodes = 104390 eff.nodes = 104390 moves = 110441 max depth = 102
_RES 1,-10.7666,-16,104390,104390,110441,102
```