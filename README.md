# DMMD (Double-matched decomposition of multi-view data)

This repository contains functions to implement DMMD method in [Dongbang Yuan and Irina Gaynanova. "Double-matched matrix decomposition for multi-view data"](https://arxiv.org/abs/2105.03396), as well as scripts to reproduce the results in the paper. 

Given two data matrices *X*<sub>1</sub> and *X*<sub>2</sub> with matched n samples, and matched p features, the method extracts additive decomposition
*X*<sub>k</sub> = *J*<sub>kc</sub> + *I*<sub>kc</sub> + *E*<sub>k</sub> = *J*<sub>kr</sub> + *I*<sub>kr</sub> + *E*<sub>k</sub>, where *J*<sub>kc</sub> is joint structure with respect to matched samples (joint column-space), *I*<sub>kc</sub> is individual structure with respect to matched samples (individual column-space), *J*<sub>kr</sub> is joint structure with respect to matched features (joint row-space), *I*<sub>kr</sub> is individual structure with respect to matched features (individual row-space), and *E*<sub>k</sub> is residual noise.

To reproduce the results, please clone this project locally and open **DMMD_Code.Rproj** in Rstudio (this makes all relative paths in scripts work as expected).

## 2. Getting started
### Main functions
*DMMDFunctions* folds contains function scripts needed for implementation of DMMD. 

**DoubleMatchedMatrixDecomposition.R** - *DMMD_v2* is the main function that performs DMMD, which includes rank selection and estimation of signal matrices given the ranks

**Double_iterative.R** - modification of DMMD to allow iterative updates of joint column and row spaces, *DMMD_i* is the main function

**FindOptMatrix.R** - *FindOpt_DM_Iterative* estimates the signal matrices **given** the ranks, this corresponds to Algorithm 1 in the paper

**Profile_Likelihood_Rank_selection.R** - contains functions for total rank estimation using profile likelihood

**Angle_Calculation.R** - contains functions for joint rank estimation using profile likelihood.

**Preliminary_Functions.R** - small supplementary functions used by the main functions above <!--. *Fnorm* function calculates Frobenius norm of a matrix. *Matscale* function does center and scale for a matrix, either row-wise or column-wise. *DoubleStandardize* function double standardizes a matrix, which makes mean 0 and variance 1 for all the rows and columns. *projection* calculates the projection matrix for a specified matrix. *GenOrthoMatrix* generates pseudo-random orthogonal matrix. *svd_recover* does rank-r svd approximation of a specified matrix.-->


### Example
To illustrate application of DMMD, we will simulate data based on DMMD model using scripts from **DoubleMatchedDataGen.R** (located in **OtherFunctions** folder)

```{r}
source("OtherFunctions/DoubleMatchedDataGen.R")
source("DMMDFunctions/Preliminary_Functions.R")
set.seed(37)
# Generate data
data = DoubleDataGen3(n = 20, p = 16, rank = c(4, 3), joint_rank_col = 2, joint_rank_row = 1, nrep = 1, std1 = 0.01, std2 = 0.01)
X1 = data$X1_list[[1]]
X2 = data$X2_list[[1]]
```

We then illustrate application of DMMD.
```{r}
# Source all necessary functions to apply DMMD
source("DMMDFunctions/Angle_Calculation.R")
source("DMMDFunctions/Profile_Likelihood_Rank_Selection.R")
source("DMMDFunctions/DoubleMatchedMatrixDecomposition.R")
source("DMMDFunctions/DMMD_iterative.R")
source("DMMDFunctions/FindOptMatrix.R")
# Apply DMMD
DMMD_result = DMMD_v2(X1, X2)
# Extract estimated ranks
DMMD_result$"Rank"
# Extract common and individual structures for 1st signal
J1c = DMMD_result$"Column Decomposition"$"Joint Column 1"
I1c = DMMD_result$"Column Decomposition"$"Individual Column 1"
J1r = DMMD_result$"Row Decomposition"$"Joint Row 1"
I1r = DMMD_result$"Row Decomposition"$"Individual Row 1"
# Verify that both decompositions give the same estimated signal
sum((J1c + I1c - J1r - I1r)^2)
# Get estimated 1st signal matrix of DMMD
A1_est = X1 - DMMD_result$Error$Error1
# Get the true signal
A1 = data$Signal1_list[[1]]
# Check the difference between estimated signal of DMMD with true signal 
sum((A1_est - A1)^2)
```
Alternatively, DMMD-i can be applied (empirically this gives slight improvement over DMMD, but has considerably higher computation cost, not recommended with large datasets). The output is a little different from DMMD as it directly returns estimated signal matrices and joint row/column spaces for convenience of later analyses.
```{r}
# Apply DMMD-i
DMMD_i_result = DMMD_i(X1, X2)
# Get estimated signal matrices of DMMD-i
A1_est_i = DMMD_i_result$A1
A2_est_i = DMMD_i_result$A2
# Get estimated basis vectors for joint subspaces
M_est = DMMD_i_result$M # joint column space (matched samples)
N_est = DMMD_i_result$N # joint row space (matched features)
# Get the joint column structure for view 1. 
J1c = M_est %*% t(M_est) %*% X1
# Get individual column structure
I1c = A1_est_i - J1c
# Get the joint row structure for view 1
J1r =  X1 %*% N_est %*% t(N_est) 
# Get individual crow structure
I1r = A1_est_i - J1r
# Check the difference between estimated signal of DMMD-i with true signal
sum((A1_est_i - A1)^2)
```

## 3. To reproduce results
DMMD algorithm does not rely on external libraries, however, for comparing different methods, corresponding libraries need to be installed. Please refer to the head of each script to make sure the packages are installed locally on your computer.

**Model_Demo_Picture** folder  - codes to generate example of DMMD decomposition on toy data from Figure 1 

**SimulationData_Setting1**, **SimulationData_Setting2**, **SimulationData_Setting3**, **Data_Setting4_FixRank**, **Data_Setting5_FixRank_SNR0.5**, **SimulationData_Setting6** folders within *Simulations*:

  - *Main.R* within each respective folder generates simulation data for corresponding setting. See Section 3 of the paper for additional details on the settings.

**RankEstimation_Setting1**, **RankEstimation_Setting2**, **RankEstimation_Setting3_Special_Final**, **Signal_Identification_Setting4**, **Signal_Identification_Setting4**, **Signal_Identification_Setting5**, **Signal_Identification_Setting6** folders within *Simulations*:

- *AJIVE.R* has implementation of [AJIVE](https://arxiv.org/pdf/1704.02060.pdf),  *SLIDE.R* has implementation of [SLIDE](https://doi.org/10.1111/biom.13108); *Main.R* runs DMMD and [JIVE](https://doi.org/10.1214/12-AOAS597), *PlotSimulations_Draft.R* combines all the results from AJIVE, SLIDE, DMMD and JIVE and plot the figures. They are saved in *Figures/Draft*.

Rank estimation comparisons also include edge distribution (ED) method in [Shu, Hai, Xiao Wang, and Hongtu Zhu. "D-CCA: A decomposition-based canonical correlation analysis for high-dimensional datasets." Journal of the American Statistical Association 115.529 (2020): 292-306](https://amstat.tandfonline.com/doi/full/10.1080/01621459.2018.1543599?casa_token=HA13MS9KztkAAAAA%3A1Q_j0Z1DWQ-32p83DDooAf1SxI318fE5HglIgRj1YyNpZY_Kv6BJ-0RTkIajA3t6vIA_QHmhuw). Also refer to [Onatski, Alexei. “Determining the number of factors from empirical distribution of eigenvalues.” The Review of Economics and Statistics, vol. 92, no. 4, 2010, pp. 1004–1016. JSTOR](https://www.jstor.org/stable/40985808?seq=1#metadata_info_tab_contents) for more information on the edge distribution method.

**New_Data_TCGA_Dim** folder  - High-dimensional simulation setting - supplement S3.4

**RankMisspecification** folder  - Signal estimation when ranks are misspecified - supplement S3.5

**Speed_test** folder - Computational comparisons -  supplement S3.6.

**Application/Soccer** folder - soccer data sets and corresponding analyses

**Application/TCGA** folder  - code for analyses of TCGA miRNA dataset (the data are in **Data** folder )

<!--References
-------
[Feng, Qing, et al. "Angle-based joint and individual variation explained." Journal of multivariate analysis 166 (2018): 241-265.](https://arxiv.org/pdf/1704.02060.pdf).

[Lock, Eric F., et al. "Joint and individual variation explained (JIVE) for integrated analysis of multiple data types." The annals of applied statistics 7.1 (2013): 523.](https://arxiv.org/pdf/1102.4110.pdf).

[Zhu, Mu, and Ali Ghodsi. "Automatic dimensionality selection from the scree plot via the use of profile likelihood." Computational Statistics & Data Analysis 51.2 (2006): 918-930.](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.90.3768&rep=rep1&type=pdf).

[Shu, Hai, Xiao Wang, and Hongtu Zhu. "D-CCA: A decomposition-based canonical correlation analysis for high-dimensional datasets." Journal of the American Statistical Association 115.529 (2020): 292-306.](https://amstat.tandfonline.com/doi/full/10.1080/01621459.2018.1543599?casa_token=HA13MS9KztkAAAAA%3A1Q_j0Z1DWQ-32p83DDooAf1SxI318fE5HglIgRj1YyNpZY_Kv6BJ-0RTkIajA3t6vIA_QHmhuw).

[Onatski, Alexei. “Determining the number of factors from empirical distribution of eigenvalues.” The Review of Economics and Statistics, vol. 92, no. 4, 2010, pp. 1004–1016. JSTOR](https://www.jstor.org/stable/40985808?seq=1#metadata_info_tab_contents).
-->
