# globalGSA - postGWAS analysis
Implementation of different Gene-set Analysis (GSA) algorithm for combining the individual pvalues of a set of genetic variants (SNPs) in a gene level p-value. 


## Simulation procedure

We generated different scenarios as described in Section 3. Data can be found in `simulatedData`. 
In order to reproduce `globalEVT`, `globalARTP` and `ARTP`, we provide a script code, `IEEE_Simulation_code.R`.
It is necessary to execute the provided code into each file folder. In addition, the execution requires the installation of the `globalGSA`, `SNPassoc`, `genetics`, `mc2d` and `multicore` packages. The new algorithm version of `globalGSA` is provided at `globalGSA_algorithm.R`.

In each scenario's folder we apply the code procedure to obtain the geneset p-value: We compute gene set p-value under `globalEVT`, `globalARTP` and `ARTP` algorithms using `globalEVT`, and `globalARTP` functions provided in the `globalGSA` package. 
File results are kept in the same directory.


Comments:
1. This code is reproduced using a multicore procedure where each independent execution is computed in a different CPU. 
You can select the total number of CPUs at the begin of the code (the default is 10). If you are unable to use multiple CPUs, this parameter must be set to "1".

2. In order to reproduce the results from the proposed simulation procedure, you must set the total number of permutations to "1,000". 


**References**

* Vilor-Tejedor N, Gonzalez JR, Calle ML. *Efficient and Powerful Method for Combining P-Values in Genome-Wide Association Studies*. IEEE/ACM Trans Comput Biol. Bioinform. 2016 Nov;13(6):1100-1106. doi: 10.1109/TCBB.2015.2509977. Epub 2015 Dec 22. PubMed PMID: 28055892.

* Vilor-Tejedor N, Calle ML. *Global adaptive rank truncated product method for gene-set analysis in association studies*. Biom J. 2014 Sep;56(5):901-11. doi:10.1002/bimj.201300192. Epub 2014 Aug 1. PubMed PMID: 25082012.

* Yu K, Li Q, Bergen AW, Pfeiffer RM, Rosenberg PS, et al. (2009). *Pathway analysis by adaptive combination of P-values*. Genetic Epidemiology 33: 700–709 [PMC free article] [PubMed]

---

Last update: 2014-11-17
Maintainer: N. Vilor-Tejedor
