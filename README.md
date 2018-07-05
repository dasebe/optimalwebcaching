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

# External libraries

This software uses the LEMON library, which is a member of the [COIN-OR](https://www.coin-or.org/) initiative.
* *Library for Efficient Modeling and Optimization in Networks*
* Copyright 2003-2012 gervary Combinatorial Optimization Research Group, EGRES
* Boost Software License, Version 1.0, see lemon/LICENSE
* Authors: lemon/AUTHORS
* Website: http://lemon.cs.elte.hu/trac/lemon
