# Derive the Optimal Caching Policy for Internet and Storage Request Traces

The tools in this repository allow calculating the optimal caching policy and its hit ratio for request traces where object sizes are variable.
More information is available in [our Sigmetrics 2018 paper](https://www.cs.cmu.edu/~dberger1/pdf/2018PracticalBound_SIGMETRICS.pdf).

### Motivation

In computer architecture, [Belady's algorithm](https://en.wikipedia.org/wiki/Page_replacement_algorithm#The_theoretically_optimal_page_replacement_algorithm) (also known as OPT or clairvoyant) gives the optimal hit ratio that can be achieved on a given trace. Clearly, this is a very useful way to benchmark the performance of caching policies.

Unfortunately, when cached objects vary in their sizes (number of bytes they take up in the cache), Belady is not anymore optimal. In fact, on real-world traces, Belady can be outperformed by recent caching policies. This is a very common case, object sizes are variable in in-memory caches like memcached, CDN caches like Varnish, and in storage systems like Ceph.

This repo introduces a new set of algorithms that enables the calculation of the optimal hit ratio and the optimal sequence of caching decisions for workloads with variable object sizes. Specifically, we relax the goal of computing OPT to obtaining accurate upper and lower bounds on OPT's hit ratio.

### Examplary Results

For variable object sizes, there are different ways of measuring a cache's performance, e.g., the object hit/miss ratio and the byte hit/miss ratio (defined below). We focus on the most common metric, the object miss ratio, and consider several online caching policies (LRU, AdaptSize, GDSF). Specifically, we compare these online policies to several bounds on the optimal cache miss ratio. We mark bounds with a U, for upper bounds, and with an L, for lower bounds.

<img alt="Object Miss Ratio of Several Online and Offline Caching Policies" src="https://raw.githubusercontent.com/dasebe/optimalwebcaching/master/OHRgoal/resultfig/optimalcachingresults.png" width="480">

This experiment shows that online caching policies can be much better than Belady. Furthermore, even advanced versions of Belady (Belady-Size) perform similarly to online policies, suggesting that online policies are near optimal with regard to their miss ratio. In contrast, our new bounds (PFOO-U and PFOO-L) show that there still remains a gap and work is now underway to bridge this gap.

# Object Hit Ratio (OHR) & Object Miss Ratio

For the object hit/miss ratio optimization goal, every object counts the same (i.e., a hit for a large 1GB object and hit for a smal 10B object will both count as a "hit").
This is appropriate in memory caches, where the cache's purpose is to minimize the number of I/O operations (random seeks) going to secondary storage.

All code for this optimization goal can be found under the directory "OHRgoal".

## Offline Algorithms

* Flow Offline Optimum (FOO): asymptotically exact derivation of optimal caching (OPT)
* Practical FOO (PFOO): fast calculation of upper and lower bounds on OPT (PFOO-U and PFOO-L)
* OFMA: prior OPT approximation from the paper "Page replacement with multi-size pages and applications to web caching"  [Irani. STOC'97]	
* LocalRatio: prior OPT approximation from the paper "A unified approach to approximating resource allocation and scheduling" [Bar-Noy, Bar-Yehuda, Freund, Naor, and Schieber. J. ACM 48 (2001)]
* various other approximations for OPT (Belady-Size, Belady, Freq-Size)

## Usage

Traces are expected in the [webcachesim](https://github.com/dasebe/webcachesim/edit/master/README.md) space-separated format with three columns (time, id, size in bytes) and a separate request on each line.

The CLI parameters of some of the tools (with examples) are as follows.

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

# Byte Hit Ratio (BHR) & Byte Miss Ratio

For the byte hit/miss ratio optimization goal, every object is weighted in proportion to the number of bytes (e.g., a hit to a 4KB object is 4x more important that a hit to a 1KB object).
This is appropriate in disk/flash caches (e.g., CDNs), where each miss incurs a bandwidth cost (which is linear in the number of missed bytes).

All code for this optimization goal will be found under the directory "BHRgoal".


# Contributors are welcome

Want to contribute? Great! We follow the [Github contribution work flow](https://help.github.com/articles/github-flow/).
This means that submissions should fork and use a Github pull requests to get merged into this code base.

This is an early-stage research prototype. There are many ways to help out.

### Bug Reports

If you come across a bug in webcachesim, please file a bug report by [creating a new issue](https://github.com/dasebe/webcachesim/issues/new). This is an early-stage project, which depends on your input!

### Write Test Cases

This project has not be thoroughly tested, test cases are likely to get a speedy merge.

### Algorithmic Issues (Network Flow) 

Both FOO and PFOO-U are much slower than they need be. See [corresponding Issue: PFOO-U is too slow](https://github.com/dasebe/optimalwebcaching/issues/1). This is fixable, but currently open.


# Academic References

We ask academic works, which built on this code, to reference the following papers:

    Practical Bounds on Optimal Caching with Variable Object Sizes
    Daniel S. Berger, Nathan Beckmann, Mor Harchol-Balter. 
    ACM SIGMETRICS, June 2018.
    Also appeared in ACM POMACS, vol. 2, issue 2, as article No. 32, June 2018.

    Towards Lightweight and Robust Machine Learning for CDN Caching
    Daniel S. Berger
    ACM HotNets, November 2018 (to appear).


# External libraries

This software uses the [LEMON](http://lemon.cs.elte.hu/trac/lemon) and [Catch2](https://github.com/catchorg/Catch2) C++ libraries.

* *LEMON: Library for Efficient Modeling and Optimization in Networks*
  * Copyright 2003-2012 gervary Combinatorial Optimization Research Group, EGRES
  * Boost Software License, Version 1.0, see lemon/LICENSE
  * Authors: lemon/AUTHORS

* *Catch2: C++ Automated Test Cases in a Header*
  * Boost Software License, Version 1.0, see tests/LICENSE.txt

