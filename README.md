# tassociation: Methods for Point and Interval Estimation of a Correlation Measure for Evaluating T-Association.

## Reference Paper

**Title:**\
Evaluation Methods for T-association of a Surrogate Endpoint

**Abstract:**\
A surrogate endpoint is a biomarker that is reasonably likely to predict clinical benefit and is used as a substitute for a direct measure of clinical benefit under the Food and Drug Administration (FDA) Accelerated Approval pathway. According to FDA guidelines, a valid surrogate endpoint must meet two associations: I-association (the association between the surrogate and true endpoints, such as disease response and overall survival) and T-association (the association between treatment effects on both endpoints, such as odds ratio and hazard ratio). I-association is commonly evaluated, but T-association is often overlooked due to the lack of appropriate statistical methods. Failure to satisfy T-association precludes a biomarker from supporting accelerated approval. To address this gap, we propose a new method to rigorously assess T-association in accordance with FDA guidelines. This method assumes that treatment effects on the surrogate and true endpoints follow a bivariate normal distribution, accounting for both within-study and between-study variances.
The key evaluation metric is the correlation coefficient, which quantifies the relationship between treatment effects on both endpoints. Model parameters, including this correlation, are estimated using maximum likelihood, restricted maximum likelihood, and a Bayesian approach. We demonstrate the method using both simulated and real-world data. The method will serve as the statistical foundation that aligns with FDA guidelines and supports future accelerated approvals. The R package to implement the proposed method is available at https://github.com/jybelindahung/T-association.

**Preprint:**\
Jo-Ying Hung, Chih-Yuan Hsu, Pei-Fang Su, and Yu Shyr (2025). Evaluation Methods for T-association of a Surrogate Endpoint. medRxiv, 2025-08. <https://www.medrxiv.org/content/10.1101/2025.08.28.25334653v2>

## Installation

### Step 1. Install and set up `cmdstanr`

``` r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)

library(cmdstanr)
```

#### Mac users

``` r
# install CmdStan
if (is.null(cmdstanr::cmdstan_path())) {
  cmdstanr::install_cmdstan()
}
```

#### Windows users
To use CmdStan on Windows, you must ensure that Rtools is installed.

(Rtools can be downloaded from: https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html)

Then run the following commands to configure the toolchain and install CmdStan:
``` r
# Check the toolchain
cmdstanr::check_cmdstan_toolchain()

# If issues are found, run automated fix
cmdstanr::check_cmdstan_toolchain(fix = TRUE)

# Re-check, should print:
# "The C++ toolchain required for CmdStan is setup properly!"
cmdstanr::check_cmdstan_toolchain()

# Then install CmdStan
cmdstanr::install_cmdstan()
```

------------------------------------------------------------------------

### Step 2. Install the `tassociation` package

``` r
install.packages("devtools")
library(devtools)
devtools::install_github("jybelindahung/T-association")
library(tassociation)
```

------------------------------------------------------------------------

## Data

To evaluate the T-association between objective response and overall survival, we collected data from six clinical trials. In each study:

-   The treatment effect on **objective response** is summarized by an **odds ratio (OR)**.\
-   The treatment effect on **overall survival** is summarized by a **hazard ratio (HR)**.\
-   The reported confidence intervals for each estimate were extracted as provided.

| Study index | Odds Ratio for Objective Response | Hazard Ratio for Overall Survival |
|:-----------:|:--------------------------------:|:--------------------------------:|
| 1           | 0.70 (95% CI 0.46-1.06)         | 1.03 (95% CI 0.86-1.24)         |
| 2           | 1.76 (95% CI 1.08-2.86)         | 0.87 (95% CI 0.71-1.05)         |
| 3           | 1.70 (95% CI 1.10-2.60)         | 0.73 (95.92% CI 0.59-0.89)      |
| 4           | 0.91 (95% CI 0.58-1.44)         | 0.76 (97.54% CI 0.56-1.02)      |
| 5           | 3.87 (95% CI 1.61-10.10)        | 0.63 (95% CI 0.42-0.93)         |
| 6           | 2.60 (95% CI 1.30-5.50)         | 0.59 (95% CI 0.44-0.79)         |

------------------------------------------------------------------------

## Usage

### Convert treatment effect estimates to meet model assumption
The function `extract_logeffect()` transforms treatment effect estimates to the log scale and computes their standard errors, which is required for the model assumptions.
This transformation is suitable for Odds Ratio and Hazard Ratio, for other treatment effect estimates, other transformations may be required.

#### Odds ratio data
```r
data.or <- data.frame(
  effect = c(0.70, 1.76, 1.70, 0.91, 3.87, 2.60),
  ci.level = rep(95, 6),
  effect.lower = c(0.46, 1.08, 1.10, 0.58, 1.61, 1.30),
  effect.upper = c(1.06, 2.86, 2.60, 1.44, 10.10, 5.50)
)

# Run the transformation
data.or.trans <- extract_logeffect(data.or)
data.or.trans
```
```
$y
[1] -0.35667494  0.56531381  0.53062825 -0.09431068  1.35325451  0.95551145

$se
[1] 0.2129625 0.2484384 0.2194431 0.2319865 0.4684528 0.3679618
```


#### Hazard ratio data
```r
data.hr <- data.frame(
  effect = c(1.03, 0.87, 0.73, 0.76, 0.63, 0.59), 
  ci.level = c(95, 95, 95.92, 97.54, 95, 95),
  effect.lower = c(0.86, 0.71, 0.59, 0.56, 0.42, 0.44), 
  effect.upper = c(1.24, 1.05, 0.89, 1.02, 0.93, 0.79)
)

# Run the transformation
data.hr.trans <- extract_logeffect(data.hr)
data.hr.trans
```
```
$y
[1] 0.0295588 -0.1392621 -0.3147107 -0.2744368 -0.4620355 -0.5276327

$se
[1] 0.09335229 0.09981828 0.10048583 0.13338983 0.20279196 0.14930331
```
The working data can be formed by:
```r
data = list(y1 = data.or.trans$y,   #log OR estimates
            se1 = data.or.trans$se, #standard errors for log OR
            y2 = data.hr.trans$y,   #log HR estimates 
            se2 = data.hr.trans$se) #standard errors for log HR
```
### Estimate by Frequentist Approach

The default frequentist method uses **Restricted Maximum Likelihood (REML)**.  
When the number of trials is fewer than 10 (e.g., 6 trials in this example), the confidence interval is automatically computed using the **parametric bootstrap**.

The `seed` argument is optional. When provided, it ensures reproducibility of the bootstrap results.

```r
t_freq(data, seed = 1127)
```
```
$Results
  Parameter   Estimates bootstrapped 95% CI
1       rho -0.79311388    (-0.980, -0.156)
2     beta1  0.40334158     (-0.081, 0.882)
3     beta2 -0.24880474    (-0.407, -0.083)
4    psi1_2  0.29543751      (0.034, 0.948)
5    psi2_2  0.02785193      (0.002, 0.097)

$Summary
[1] "Point estimates obtained using Restricted Maximum Likelihood (REML). 
     95% confidence intervals computed using the bootstrap method with 1000 bootstrap replicates."
```

The estimated correlation (rho) in this example using the frequentist approach is −0.79, with a 95% bootstrap confidence interval of (−0.980, −0.156). 

To specify the point estimation method (`ML` or `REML`) and the interval estimation method (`normal` or `bootstrap`), use:

```r
t_freq(data.oros, method = "ML",   interval.method = "normal")
t_freq(data.oros, method = "ML",   interval.method = "bootstrap")
t_freq(data.oros, method = "REML", interval.method = "normal")
```
### Estimate by Bayesian approach
This `t_bayes` function estimates the T-association using a Bayesian approach with MCMC, returning posterior summaries for all parameters.
The default setting uses a **half-$t(3, 0.3)$** prior on $\psi$ and a **Normal(0, 1)** prior on the Fisher-$z$ transformed $\rho$.

The `seed` argument is optional. When provided, it ensures reproducibility of the MCMC results.

```r
t_bayes(data, seed = 1127)
```
```
$Results
  Parameter Posterior Mean Posterior SD          95% CrI
1       rho         -0.556        0.261  (-0.904, 0.057)
2     beta1          0.378        0.251   (-0.09, 0.929)
3     beta2         -0.243        0.098 (-0.444, -0.055)
4    psi1_2          0.275        0.272   (0.028, 0.959)
5    psi2_2          0.040        0.050   (0.001, 0.165)


$prior_summary
[1] "Psi prior: half-t(3, 0.3) on SD (psi)"
[2] "Rho prior: fisherz"
```

The estimated correlation (rho) in this example using the Bayesian approach is -0.556, with a 95% credible interval of (-0.904, 0.057). 

### Adjusting Priors

You may adjust the priors used for the variance parameters (`psi`) and the correlation parameter (`rho`).  
The available options and default values are summarized below.

#### Correlation parameter (`rho`)

Two prior choices are available:

- **Normal(0, 1)** prior on the Fisher-z transformed `rho`
- **Uniform(-1, 1)** prior directly on `rho`

#### Variance parameters (`psi`)

Two prior families are supported, each with customizable hyperparameters:

| Prior type                  | Default parameters               | Customize via                |
|-----------------------------|----------------------------------|------------------------------|
| Half-t on `psi`             | df = 3, scale = 0.3              | `half_t_df`, `half_t_scale` |
| Inverse-Gamma on `psi^2`    | shape = 0.001, scale = 0.001     | `ig_shape`, `ig_scale`      |

#### Examples, 
```r
# Default half-t prior: half-t(3,0.3) on psi; Uniform(-1,1) prior on rho
t_bayes(data, psi.prior.dist = "half-t", rho.prior.dist = "uniform")

# Default inverse-gamma prior: IG(0.001, 0.001) on psi^2; N(0,1) prior on Fisher-z transformed rho
t_bayes(data, psi.prior.dist = "inverse-gamma", rho.prior.dist = "fisherz")

# Half-t prior with custom parameters; Uniform(-1,1) prior on rho
t_bayes(
  data,
  psi.prior.dist = "half-t",
  half_t_df = 3,
  half_t_scale = 0.5,
  rho.prior.dist = "uniform"
)
```

------------------------------------------------------------------------
