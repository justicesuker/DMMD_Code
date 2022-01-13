# DMMD_Code

Description 
-------
Code for reproducing results shown in "Double-matched matrix decomposition for multi-view data". Please run **DMMD_Code.Rproj** before opening individual files to make sure the working directory is correct.

## Main functions
Listed below are functions that are used to perform double-matched matrix decomposition. They are all located in *Simulations/MyFunction*.

**Preliminary_Functions.R** - Some small convenient preliminary functions. *Fnorm* function calculates Frobenius norm of a matrix. *Matscale* function does center and scale for a matrix, either row-wise or column-wise. *DoubleStandardize* function double standardizes a matrix, which makes mean 0 and variance 1 for all the rows and columns. *projection* calculates the projection matrix for a specified matrix. *GenOrthoMatrix* generates pseudo-random orthogonal matrix. *svd_recover* does rank-r svd approximation of a specified matrix.

**Profile_Likelihood_Rank_selection.R** - The main function is *ProfileLikCluster*, which separates an input vector into two groups using profile likelihood method. It could be used to estimate rank by giving the singular values of the matrix.

**Select_ED_Rank.R** - *Select_ED_Rank* function uses edge distribution (ED) method to calculate the rank i.e., the number of factors. The function is translated from the python code *Select_ED_Rank_K* in D-CCA.py from the supplementary material of [Shu, Hai, Xiao Wang, and Hongtu Zhu. "D-CCA: A decomposition-based canonical correlation analysis for high-dimensional datasets." Journal of the American Statistical Association 115.529 (2020): 292-306](https://amstat.tandfonline.com/doi/full/10.1080/01621459.2018.1543599?casa_token=HA13MS9KztkAAAAA%3A1Q_j0Z1DWQ-32p83DDooAf1SxI318fE5HglIgRj1YyNpZY_Kv6BJ-0RTkIajA3t6vIA_QHmhuw). Also refer to [Onatski, Alexei. “Determining the number of factors from empirical distribution of eigenvalues.” The Review of Economics and Statistics, vol. 92, no. 4, 2010, pp. 1004–1016. JSTOR](https://www.jstor.org/stable/40985808?seq=1#metadata_info_tab_contents) for more information on the edge distribution method.

**Angle_Calculation.R** - There are two functions in the file. One is *angle_cal*, which calculates principal angles between column spaces of two matrices. The other is *joint_angle_cluster* that estimates joint rank by method of profile likelihood.

**FindOptMatrix.R** - *FindOpt_DM_Iterative* solves the optimization problem of (1) written in the manuscript.

**DoubleMatchedMatrixDecomposition.R** - *DMMD_v2* is the main function of the project, which performs double-matched matrix decomposition.

**DoubleMatchedDataGen.R** - *DoubleDataGen3* is the function that generates double matched data used for simulation.

## Simulation data
Six settings of simulation data are generated for comparing the performance of DMMD with other competitors on rank estimation as well as signal identification. Refer to Section 3 in the draft for more details. They are located in *Simulations*.

**SimulationData_Setting1** - *Main.R* generates the simulation data for Setting 1 and save the data in *Data.RData*.

**SimulationData_Setting2** - *Main.R* generates the simulation data for Setting 2 and save the data in *Data.RData*.

**SimulationData_Setting3** - *Main.R* generates the simulation data for Setting 3 and save the data in *Data.RData*.

**Data_Setting4_FixRank** - *Main.R* generates the simulation data for Setting 4 and save the data in *Data.RData*.

**Data_Setting5_FixRank_SNR0.5** - *Main.R* generates the simulation data for Setting 5 and save the data in *Data.RData*.

**SimulationData_Setting6** - *Main.R* generates the simulation data for Setting 6 and save the data in *Data1.RData* and *Data2.RData*.

## Simulations scripts
Different methods are run on the simulated data on six settings. The results and figures are saved. They are located in *Simulations*.

**RankEstimation_Setting1** - *AJIVE.R* runs AJIVE on simulation data for Setting 1 and saves the result in *AJIVEoutput.RData*; *SLIDE.R* runs SLIDE on simulation data for Setting 1 and saves the result in *SLIDEoutput.RData*; *Main.R* runs DMMD and JIVE on simulation data for Setting 1 and saves the result in *output.RData*, *PlotSimulations_Draft.R* combines all the results from AJIVE, SLIDE, DMMD and JIVE and plot the figures. They are saved in *Figures/Draft*.

**RankEstimation_Setting2** - The files have Similar functions as **RankEstimation_Setting1**. The results are based on simulation data for Setting 2.

**RankEstimation_Setting3_Special_Final** - The files have Similar functions as **RankEstimation_Setting1**. The results are based on simulation data for Setting 3.

**Signal_Identification_Setting4** - The files have Similar functions as **RankEstimation_Setting1**. The results are based on simulation data for Setting 4.

**Signal_Identification_Setting5** - The files have Similar functions as **RankEstimation_Setting1**. The results are based on simulation data for Setting 5.

**Signal_Identification_Setting6_Special** - The files have Similar functions as **RankEstimation_Setting1**. The results are based on simulation data for Setting 6.

## Demonstration figure 
The **Model_Demo_Picture** folder contains the code expressing how Figure 1 in the draft is generated.

## Application
**Soccer** folder contains the soccer data sets we use for analysis and the analysis results shown in draft. 
**TCGA** folder contains the miRNA data sets we use for analysis and the analysis results shown in draft. 

## Example
```{r}
source("Simulations/MyFunction/Angle_Calculation.R")
source("Simulations/MyFunction/Profile_Likelihood_Rank_Selection.R")
source("Simulations/MyFunction/DoubleMatchedMatrixDecomposition.R")
source("Simulations/MyFunction/FindOptMatrix.R")
source("Simulations/MyFunction/Preliminary_Functions.R")
source("Simulations/MyFunction/Select_ED_Rank.R")
source("Simulations/MyFunction/DoubleMatchedDataGen.R")
source("Simulations/MyFunction/DMMD_Irina.R")
set.seed(37)
data = DoubleDataGen3(n = 20, p = 16, rank = c(5, 4), joint_rank_col = 2, joint_rank_row = 1, nrep = 1, std1 = 0.01, std2 = 0.01)
X1 = data$X1_list[[1]]
X2 = data$X2_list[[1]]
DMMD_result = DMMD_v2(X1, X2)
DMMD_i_result = DMMD_Irina_unknown(X1, X2)
DMMD_i_result2 = DMMD_Irina(X1,X2, r1 = 8, r2 = 6, rc = 1, rr = 1)
A1_est = X1 - DMMD_result$Error$Error1
A1 = data$Signal1_list[[1]]
head(A1_est - A1)
```

References
-------
[Feng, Qing, et al. "Angle-based joint and individual variation explained." Journal of multivariate analysis 166 (2018): 241-265.](https://arxiv.org/pdf/1704.02060.pdf).

[Lock, Eric F., et al. "Joint and individual variation explained (JIVE) for integrated analysis of multiple data types." The annals of applied statistics 7.1 (2013): 523.](https://arxiv.org/pdf/1102.4110.pdf).

[Zhu, Mu, and Ali Ghodsi. "Automatic dimensionality selection from the scree plot via the use of profile likelihood." Computational Statistics & Data Analysis 51.2 (2006): 918-930.](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.90.3768&rep=rep1&type=pdf).

[Shu, Hai, Xiao Wang, and Hongtu Zhu. "D-CCA: A decomposition-based canonical correlation analysis for high-dimensional datasets." Journal of the American Statistical Association 115.529 (2020): 292-306.](https://amstat.tandfonline.com/doi/full/10.1080/01621459.2018.1543599?casa_token=HA13MS9KztkAAAAA%3A1Q_j0Z1DWQ-32p83DDooAf1SxI318fE5HglIgRj1YyNpZY_Kv6BJ-0RTkIajA3t6vIA_QHmhuw).

[Onatski, Alexei. “Determining the number of factors from empirical distribution of eigenvalues.” The Review of Economics and Statistics, vol. 92, no. 4, 2010, pp. 1004–1016. JSTOR](https://www.jstor.org/stable/40985808?seq=1#metadata_info_tab_contents).


