# HMDBN
Hidden Markov induced Dynamic Bayesian Network for recovering time evolving gene regulatory networks

## Description
Dynamic Bayesian Networks (DBN) have been widely used to recover gene regulatory relationships from time-series data in computational systems biology. Its standard assumption is `stationarity`, and therefore, several research efforts have been recently proposed to relax this restriction. However, those methods suffer from three challenges: `long running time`, `low accuracy` and `reliance on parameter settings`. To address these problems, we propose a novel `non-stationary DBN model` by extending each hidden node of Hidden Markov Model into a DBN (called `HMDBN`), which properly handles the underlying time-evolving networks. It has the following attractive contributions:
-    An improved structural EM algorithm is proposed to learn the HMDBN. It dramatically reduces searching space, thereby substantially improving computational efficiency. 
-    A novel generalized Bayesian Information Criterion under the non-stationary assumption (called `BWBIC`) is derived, which can help significantly improve the reconstruction accuracy and largely reduce over-fitting. 
-    The re-estimation formulas for all parameters of our model are derived, enabling us to avoid reliance on parameter settings. 
-    Compared to the state-of-the-art methods, the experimental evaluation of our proposed method on both synthetic and real biological data demonstrates more stably high prediction accuracy and significantly improved computation efficiency, even with no prior knowledge and parameter settings.

## Depends
- Matlab

## Tutorial
See [HMDBN tutorial](https://github.com/zhushijia/HMDBN/blob/master/man/HMDBN_MatlabDoc.pdf)

## Citation
Shijia Zhu & Yadong Wang, Hidden Markov induced Dynamic Bayesian Network for recovering time evolving gene regulatory networks ,Scientific Reports 5: 17841 (2015). [(link)](https://www.nature.com/articles/srep17841.pdf)
