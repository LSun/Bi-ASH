# Bi-ASH: Covariate-modulated Adaptive Shrinkage by Empirical Bayes

We have observations $\hat\beta_j$ of true effect $\beta_j$ with estimated measurement error $\hat s_j$, as well as covariate $X_j$.  The true effects presumably have stochastic orders according to their covariates.  We hope to make covariate-modulated inference of the true effects by sharing information among all the data using empirical Bayes methods.

$$
\begin{array}{c}
\beta_j = \alpha_j + X_j\gamma_j\\
\hat\beta_j|\hat s_j \sim N(\beta_j,\hat s_j^2) \text{ independently}\\
j=1,\ldots,J
\end{array}
$$

where

$$
\begin{array}{c}
\alpha_j \sim f^{\alpha}(\cdot) \text{ exchangeably, }f^\alpha \text{ uniform, maybe symmetric}\\
\gamma_j \sim f^{\gamma}(\cdot) \text{ exchangeably, }f^\gamma \text{ uniform, maybe symmetric}\\
X_j \text{ binary or continuous}
\end{array}
$$

In this setting, if $X_j$ is continuous, larger $|X_j|$ leads to stochastically larger $\beta_j$; if $X_j$ is binary (0 or 1), true effects in Group1 is stochastically larger than those in Group0.