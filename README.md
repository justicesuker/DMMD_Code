# DMMD_Code

Code for reproducing results and corresponding graphs shown in "Double-matched matrix decomposition for multi-view data". 

Description 
-------

We consider the case of two matched datasets, where the matching is performed both by rows and by columns. This situation, arises, for example, when the expression levels of the same genes are measured across two tissues of the same subject. The goal of the analysis is to identify common as well as individual signals within each tissue, where common signals can be separated into two parts: common across subjects, and common across genes. In this work, we combine the view of row and column directions, then prove a lemma that guarantees the existence and uniqueness of matrix factorization in double matched cases. Based on the lemma, a new technique which could calculate the common and individual signals in both ways is introduced. In the section of algorithm, we provide a valid way to calculate the decomposition efficiently. After that, we provide the theoretical guarantee of convergence on the algorithm. In the end, we compare our methods with other methods on performance using simulated data as well as real data.

## Main functions
Listed below are functions that are used to perform double-matched matrix decomposition. They are all located in *Simulations/MyFunction*

**Preliminary_Functions.R** - Some small convenient preliminary functions. *Fnorm* function calculates Frobenius norm of a matrix. *Matscale* function does center and scale for a matrix, either row-wise or column-wise. *DoubleStandardize* function double standardizes a matrix, which makes mean 0 and variance 1 for all the rows and columns. *projection* calculates the projection matrix for a specified matrix. *GenOrthoMatrix* generates pseudo-random orthogonal matrix. *svd_recover* does rank-r svd approximation of a specified matrix.

**Profile_Likelihood_Rank_selection.R** - The main function is *ProfileLikCluster*, which separates an input vector into two groups using profile likelihood method. It could be used to estimate rank by giving the singular values of the matrix.

**Select_ED_Rank.R** - *Select_ED_Rank* function uses edge distribution (ED) method to calculate the rank i.e., the number of factors. The function is translated from the python code *Select_ED_Rank_K* in D-CCA.py from the supplementary material of [Shu, Hai, Xiao Wang, and Hongtu Zhu. "D-CCA: A decomposition-based canonical correlation analysis for high-dimensional datasets." Journal of the American Statistical Association 115.529 (2020): 292-306](https://amstat.tandfonline.com/doi/full/10.1080/01621459.2018.1543599?casa_token=HA13MS9KztkAAAAA%3A1Q_j0Z1DWQ-32p83DDooAf1SxI318fE5HglIgRj1YyNpZY_Kv6BJ-0RTkIajA3t6vIA_QHmhuw). Also refer to [Onatski, Alexei. “Determining the number of factors from empirical distribution of eigenvalues.” The Review of Economics and Statistics, vol. 92, no. 4, 2010, pp. 1004–1016. JSTOR](https://www.jstor.org/stable/40985808?seq=1#metadata_info_tab_contents) for more information on the edge distribution method.

**Angle_Calculation.R** - There are two functions in the file. One is *angle_cal*, which calculates principal angles between column spaces of two matrices. The other is *joint_angle_cluster* that estimates joint rank by method of profile likelihood.

**FindOptMatrix.R** - *FindOpt_DM_Iterative* solves the optimization problem of (1) written in the manuscript.

**DoubleMatchedMatrixDecomposition.R** - *DMMD_v2* is the main function of the project, which performs double-matched matrix decomposition.

**DoubleMatchedDataGen.R** - *DoubleDataGen2* is the function that generates double matched data used for simulation.

**DMMD_Impute.R** - The main function is *DMMD_Impute* which performs an iterative procedure to do imputation based on DMMD method.

## Simulation data
Three settings (currently two, the potential third one will be added soon) of data are generated for comparing DMMD method with other competitors on rank estimation as well as signal identification. Refer to section 3 in the draft manuscript. Folders below are located in the folder of *Simulations*

**SimulationData** - The folder contains the data for simulation setting 1 mentioned in the draft manuscript. *Main.R* is the program that generates the data. *Data.RData* stores the generated data.

**SpecialRankData** - The folder contains the data for simulation setting 2 mentioned in the draft manuscript. *Main.R* is the program that generates the data. *Data.RData* stores the generated data.

**TBD** - To be added. Signal to noise ratio is 0.5.

## Simulations scripts

**RankEstimation_With_ED** - Simulation setting 1, compare the performance on rank estimation. 

**Signal_Relative_20200604** - Simulation setting 1, compare the performance on signal identification.

**SpecialRank_20200604_ED** - Simulation setting 2, compare the performance on rank estimation. 

**SpecialRank_RelativeSignal_20200604** - Simulation setting 2, compare the performance on signal identification.

**TBD_Rank_Setting3** - Simulation setting 3, compare rank estimation performance.

**TBD_Signal_Setting3** - Simulation setting 3, compare the performance on signal identification.

References
-------
[Feng, Qing, et al. "Angle-based joint and individual variation explained." Journal of multivariate analysis 166 (2018): 241-265.](https://arxiv.org/pdf/1704.02060.pdf).

[Lock, Eric F., et al. "Joint and individual variation explained (JIVE) for integrated analysis of multiple data types." The annals of applied statistics 7.1 (2013): 523.](https://arxiv.org/pdf/1102.4110.pdf).

[Zhu, Mu, and Ali Ghodsi. "Automatic dimensionality selection from the scree plot via the use of profile likelihood." Computational Statistics & Data Analysis 51.2 (2006): 918-930.](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.90.3768&rep=rep1&type=pdf).

[Shu, Hai, Xiao Wang, and Hongtu Zhu. "D-CCA: A decomposition-based canonical correlation analysis for high-dimensional datasets." Journal of the American Statistical Association 115.529 (2020): 292-306.](https://amstat.tandfonline.com/doi/full/10.1080/01621459.2018.1543599?casa_token=HA13MS9KztkAAAAA%3A1Q_j0Z1DWQ-32p83DDooAf1SxI318fE5HglIgRj1YyNpZY_Kv6BJ-0RTkIajA3t6vIA_QHmhuw).

[Onatski, Alexei. “Determining the number of factors from empirical distribution of eigenvalues.” The Review of Economics and Statistics, vol. 92, no. 4, 2010, pp. 1004–1016. JSTOR](https://www.jstor.org/stable/40985808?seq=1#metadata_info_tab_contents).

