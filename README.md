# Derive the Optimal Caching Policy for Internet and Storage Request Traces

The tools in this repository allow calculating the optimal caching policy and its hit ratio for request traces where object sizes are variable.
More information is available in [our Sigmetrics 2018 paper](https://www.cs.cmu.edu/~dberger1/pdf/2018PracticalBound_SIGMETRICS.pdf).

    Practical Bounds on Optimal Caching with Variable Object Sizes
    Daniel S. Berger, Nathan Beckmann, Mor Harchol-Balter. 
    ACM SIGMETRICS, June 2018.
    Also appeared in ACM POMACS, vol. 2, issue 2, as article No. 32, June 2018.

## Overview

* Flow Offline Optimum (FOO): asymptotically exact derivation of optimal caching (OPT)
* Practical FOO (PFOO): fast calculation of upper and lower bounds on OPT (PFOO-U and PFOO-L)
* OFMA: prior OPT approximation from the paper "Page replacement with multi-size pages and applications to web caching"  [Irani. STOC'97]	
* LocalRatio: prior OPT approximation from the paper "A unified approach to approximating resource allocation and scheduling" [Bar-Noy, Bar-Yehuda, Freund, Naor, and Schieber. J. ACM 48 (2001)]
* various other approximations for OPT (Belady-Size, Belady, Freq-Size)

## Usage

Traces are expected in the [webcachesim](https://github.com/dasebe/webcachesim/edit/master/README.md) space-separated format with three columns (time, id, size in bytes) and a separate request on each line.

The CLI parameters of the tools (with examples) are as follows.

* FOO 
  * format (four parameters, all required) and example:
    ```
    ./foo [trace name] [cachesize in bytes] [pivot rule] [output name]
    ./foo trace.txt 1073741824 4 foo_decision_variables.txt
    ``` 
  * pivot rule denotes the [network simplex's pivot rule](http://lemon.cs.elte.hu/pub/doc/latest/a00269.html)

* PFOO-U
  * format (four parameters, all required) and example:
    ```
    ./pfoou [trace name] [cachesize in bytes] [pivot rule] [step size] [output name]
    ./pfoou trace.txt 1073741824 4 500000 pfoo_decision_variables.txt
    ``` 
  * step size denotes the number of decision variables that PFOO-U changes in each iteration, 500k is a good starting point. (Lower is faster, higher has better accuracy)
  
* PFOO-L
  * format (four parameters, all required) and example:
    ```
    ./pfool [trace name] [cachesize in bytes] [output name]
    ./pfoou trace.txt 1073741824 pfoo_decision_variables.txt
    ``` 


## Contributors are welcome

Want to contribute? Great! We follow the [Github contribution work flow](https://help.github.com/articles/github-flow/).
This means that submissions should fork and use a Github pull requests to get merged into this code base.

This is an early-stage research prototype. There are many ways to help out.

### Bug Reports

If you come across a bug in webcachesim, please file a bug report by [creating a new issue](https://github.com/dasebe/webcachesim/issues/new). This is an early-stage project, which depends on your input!

### Write test cases

This project has not be thoroughly tested, test cases are likely to get a speedy merge.

### Algorithmic contributions

Both FOO and PFOO-U are much slower than they need be. See [corresponding Issue: PFOO-U is too slow](https://github.com/dasebe/optimalwebcaching/issues/1). This is fixable, but I do not have time to work on this right now.

# External libraries

This software uses the [LEMON](http://lemon.cs.elte.hu/trac/lemon) and [Catch2](https://github.com/catchorg/Catch2) C++ libraries.

* *LEMON: Library for Efficient Modeling and Optimization in Networks*
  * Copyright 2003-2012 gervary Combinatorial Optimization Research Group, EGRES
  * Boost Software License, Version 1.0, see lemon/LICENSE
  * Authors: lemon/AUTHORS

* *Catch2: C++ Automated Test Cases in a Header*
  * Boost Software License, Version 1.0, see tests/LICENSE.txt

