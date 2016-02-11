---
author: Bryce Bartlett
date: 2/7/2016
title: "Outline of Hierarchical Bayesian Changepoint Model."
tags: Bayesian, time series
---

#Overview#

This outlines the equation and provides an explanaiton of why I'm using a modified Bayeslian Changepoint model.

#Summary#

A classic explanatino of a Bayesian changepoint model was described by Carlin, Gefland, and Smith in 1992 ("Hierarchical Bayesian Analysis of Cahgnepoint Problems"). This sort of model is a bugs example model. The basic concept idea of the changepointmodel is that given a time series, there are finite nubmer of distributions that explain the data generation process, and these distributions change at various points in time (p. 390-391). Carlin, Gefland, and Smith outline a generalizable linear model to identify the changepoint between log-linear rates for the flow of water. The model I use is fundamentally the same, but with a *known* changepoint, starting in 1999, and includes a hierarchcial structure for the particular cells used to calculate the rates. The basic model is as follows:

Level 1, Within Cell (Likelihood):

**Check Lynch on structure; think about whether the (c) should go below and the 1, 2 should go above...**
**Can you specify priors as equal like you have in the model ??**

$$
y_{ct} \sim \begin{cases}
   N(\alpha^{(c)}_1 + \beta_1 t, \sigma^2_1) \mbox{, if } t<0  \\
   N(\alpha^{(c)}_2 + \beta_2 t, \sigma^2_2) \mbox{, if } t \geq 0   
\end{cases}
$$

Level 2, Between Cell (Prior):

$$
\begin{align}
\begin{bmatrix}
  \alpha_1 \\  \alpha_2 
\end{bmatrix}
\sim \text{MVN} \Big(
\begin{bmatrix}
 \mathbf{Z_c \Gamma_1} \\  \mathbf{Z_c \Gamma_2}
\end{bmatrix}
,\mathbf{\Sigma} \Big)
\end{align}
$$

Keep the above in the model; put the below in a technical supplement.

$$
\beta \sim N(0,5)
$$

**note that LKJ is described on pp. 582 and 576 of BDAIII. And p. 437 includes a comment on the half-cauchy**     

Hyperprior:

$$
\begin{aligned}
& \gamma \sim \text{N}(0,5) \\ 
& \sigma \sim N(0,5) \\
& \mathbf{\Sigma} = Diag(\xi) \mathbf{\Omega} Diag(\xi) \\
& \text{where  } 
\begin{align}
& \mathbf{\Omega} \sim \text{LkjCorr}(1)  \\
 & \xi \sim \text{Half-Cauchy}(5)
 \end{align}
\end{aligned}
$$


key points: (1) I have a known changepoint; (2) takes care of heteroskedasticity; (3) takes care of serial errors; and (4) hierarchical is better as fully bayesian. Bayesian also gives additional information. Fundamentally similar to a dummy interaction across the whole model with robust s.e. in frequentist (add one in if you can in an appendix or say it is available on request).