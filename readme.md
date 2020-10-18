Overview
--------

The **REsmcure** package provides an extension to the CRAN-package
[smcure](https://github.com/cran/smcure) to fit Semi-Parametric (Cox
Proportional Hazards) Mixture Cure models with Bivariate Gaussian random
effects (one both the incidence and the latency part of the model). The
code allows for (and estimates) correlation of both random effects. The
model estimation is based on a paper by Lai & Yau (2008).

Content of the code
-------------------

-   simsomdata is a function that can be used to simulate some survival
    data with (correlated) random effects.

-   REsmsurv extends the smsurv function from the **smcure** package so
    it takes the presence of random effects into account. This function
    computes the baseline survival function for given data and (initial)
    parameter values. *for internal use*

-   REem implements the Expectation-Maximization (with Inner and Outer
    loop) for parameter estimation. The outer loop updates the fixed
    parameters, the BLUPs for the random effects and the baseline
    survival function, while the inner loop updates the parameters of
    the random effects’ distribution. *for internal use*

-   REsmcure is the wrapper function to estimate a mixture cure survival
    model with correlated random effects. The function allows different
    fixed effect models on the incidence (*cureform*) and the latency
    (*formula*). The random effect identifier is named (variable name as
    a string) using the RE argument.

-   sebootjack is a post-hoc function to add robust standard errors to
    the model object. Robust standard errors are obtained using either a
    parametric bootstrap procedure or a leave-k-out jackknife approach
    on the given dataset.

-   print.REsmcure is the printing function to obtain formatted model
    output.

How to use the code
-------------------

``` r
#devtools::install_github('HansTierens/MMixCure',dep=T)
library(MMixCure)
```

### Simulate some data

``` r
i=123456789
n=1000
ng=10
coefficients<-list(gamma0 = c(0,.40,-.30,-.30,.25),  beta0 = c(0,.45,.35,-.25,-.15))
revcovmat<-list(varcure=1,  varfrail=1,  corr=-.8)

data<-simsomdata(i,n,ng,coefficients,revcovmat)
knitr::kable(data[1:10,])
```

|          Y|  Delta|          Xs|          Xd|   Ds|   Dd|  grp|
|----------:|------:|-----------:|-----------:|----:|----:|----:|
|  6.1558647|      0|  -0.7339469|  -1.3956548|    0|    0|    5|
|  1.0461923|      1|   1.3057519|  -1.6428802|    0|    0|    6|
|  7.0937835|      0|  -1.3405564|  -0.8882187|    0|    1|    9|
|  1.7904848|      1|  -0.2297080|   0.6969024|    0|    1|    9|
|  2.2355639|      0|   0.3829667|   1.5881038|    1|    0|    6|
|  0.3325948|      1|   2.5447458|   1.3653490|    1|    0|    4|
|  0.4227490|      1|  -0.9422011|  -0.8880386|    0|    0|    1|
|  0.9253535|      1|  -0.7104432|   0.4318116|    0|    0|    8|
|  0.6459655|      0|   0.4594213|  -0.4245565|    1|    1|   10|
|  2.4323292|      0|   0.0190030|   0.1624558|    0|    0|    6|

``` r
KMest<-survfit(coxph(Surv(Y,Delta)~1,data))
plot(KMest,main='Simulated Survival Data',ylab='Survival Prob.',xlab='Time')
```

![](readme_files/figure-markdown_github/simdata-1.png)

### Fit the Model

``` r
coxfit<-coxph(Surv(Y,Delta)~1+Xs+Xd+Ds+Dd,data=data,ties='breslow')
summary(coxfit)
```

    ## Call:
    ## coxph(formula = Surv(Y, Delta) ~ 1 + Xs + Xd + Ds + Dd, data = data, 
    ##     ties = "breslow")
    ## 
    ##   n= 1000, number of events= 458 
    ## 
    ##        coef exp(coef) se(coef)      z Pr(>|z|)    
    ## Xs  0.41191   1.50970  0.04774  8.628  < 2e-16 ***
    ## Xd -0.05433   0.94711  0.04793 -1.134  0.25691    
    ## Ds -0.25255   0.77682  0.09368 -2.696  0.00702 ** 
    ## Dd  0.07760   1.08069  0.09361  0.829  0.40714    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##    exp(coef) exp(-coef) lower .95 upper .95
    ## Xs    1.5097     0.6624    1.3748    1.6578
    ## Xd    0.9471     1.0558    0.8622    1.0404
    ## Ds    0.7768     1.2873    0.6465    0.9334
    ## Dd    1.0807     0.9253    0.8995    1.2983
    ## 
    ## Concordance= 0.636  (se = 0.013 )
    ## Likelihood ratio test= 82.24  on 4 df,   p=<2e-16
    ## Wald test            = 81.99  on 4 df,   p=<2e-16
    ## Score (logrank) test = 82.66  on 4 df,   p=<2e-16

``` r
frailtyfit<-coxph(Surv(Y,Delta)~1+Xs+Xd+Ds+Dd+frailty.gaussian(grp),data=data,ties='breslow')
summary(frailtyfit)
```

    ## Call:
    ## coxph(formula = Surv(Y, Delta) ~ 1 + Xs + Xd + Ds + Dd + frailty.gaussian(grp), 
    ##     data = data, ties = "breslow")
    ## 
    ##   n= 1000, number of events= 458 
    ## 
    ##                       coef     se(coef) se2     Chisq DF   p      
    ## Xs                     0.40741 0.04767  0.04765 73.05 1.00 1.3e-17
    ## Xd                    -0.05599 0.04847  0.04843  1.33 1.00 2.5e-01
    ## Ds                    -0.25268 0.09486  0.09470  7.10 1.00 7.7e-03
    ## Dd                     0.09481 0.09398  0.09394  1.02 1.00 3.1e-01
    ## frailty.gaussian(grp)                           54.62 8.45 8.3e-09
    ## 
    ##    exp(coef) exp(-coef) lower .95 upper .95
    ## Xs    1.5029     0.6654    1.3689    1.6501
    ## Xd    0.9455     1.0576    0.8599    1.0398
    ## Ds    0.7767     1.2875    0.6449    0.9354
    ## Dd    1.0994     0.9096    0.9145    1.3218
    ## 
    ## Iterations: 6 outer, 23 Newton-Raphson
    ##      Variance of random effect= 0.1515622 
    ## Degrees of freedom for terms= 1.0 1.0 1.0 1.0 8.4 
    ## Concordance= 0.66  (se = 0.66 )
    ## Likelihood ratio test= 149.6  on 12.44 df,   p=<2e-16

``` r
spmcfit<-smcure(Surv(Y,Delta)~1+Xs+Xd+Ds+Dd,cureform=~1+Xs+Xd+Ds+Dd,data=data,model='ph',link='logit',Var=TRUE)
```

    ## Program is running..be patient... done.
    ## Call:
    ## smcure(formula = Surv(Y, Delta) ~ 1 + Xs + Xd + Ds + Dd, cureform = ~1 + 
    ##     Xs + Xd + Ds + Dd, data = data, model = "ph", link = "logit", 
    ##     Var = TRUE)
    ## 
    ## Cure probability model:
    ##                Estimate  Std.Error    Z value     Pr(>|Z|)
    ## (Intercept)  0.33373267 0.13125134  2.5426991 1.099999e-02
    ## Xs           0.44518469 0.09195094  4.8415457 1.288331e-06
    ## Xd          -0.19522292 0.08423294 -2.3176552 2.046807e-02
    ## Ds          -0.24432276 0.17725276 -1.3783862 1.680841e-01
    ## Dd           0.02803339 0.18260473  0.1535195 8.779886e-01
    ## 
    ## 
    ## Failure time distribution model:
    ##      Estimate  Std.Error    Z value     Pr(>|Z|)
    ## Xs  0.3402169 0.07540967  4.5115816 6.434603e-06
    ## Xd  0.1541455 0.05166311  2.9836667 2.848168e-03
    ## Ds -0.2720835 0.13331783 -2.0408637 4.126438e-02
    ## Dd  0.0663760 0.12442585  0.5334583 5.937164e-01

``` r
fit<-REsmcure(Surv(Y,Delta)~1+Xs+Xd+Ds+Dd,cureform=~1+Xs+Xd+Ds+Dd,RE='grp',rho=T, data=data, emmax = 50, eps = 1e-07) 
print.REsmcure(fit)
```

    ## 
    ##  Call:
    ## REsmcure(formula = Surv(Y, Delta) ~ 1 + Xs + Xd + Ds + Dd, cureform = ~1 + 
    ##     Xs + Xd + Ds + Dd, RE = "grp", rho = T, data = data, emmax = 50, 
    ##     eps = 1e-07)
    ## 
    ## Cure probability model:
    ##             Estimate Std.Error  Z value Pr(>|Z|)
    ## (Intercept)  0.03318   0.33667  0.09855  0.92150
    ## Xs           0.45936   0.07305  6.28847  0.00000
    ## Xd          -0.29312   0.07474 -3.92195  0.00009
    ## Ds          -0.42290   0.14392 -2.93835  0.00330
    ## Dd           0.04380   0.14240  0.30762  0.75837
    ## 
    ## 
    ## Failure time distribution model:
    ##    Estimate Std.Error  Z value Pr(>|Z|)
    ## Xs  0.44312   0.05034  8.80326  0.00000
    ## Xd  0.32353   0.04972  6.50725  0.00000
    ## Ds -0.15921   0.09814 -1.62218  0.10477
    ## Dd -0.00897   0.09571 -0.09370  0.92535
    ## 
    ## Random Effects:
    ##                  Variance Correlation
    ## Cure/Incidence       0.98            
    ## Survival/Latency     1.33      -0.841

Citations and References
------------------------

-   This package builds on the **smcure** package, so please do refer
    to:

Cai, C., Zou, Y., Peng, Y., & Zhang, J. (2012). smcure: An R-package for
estimating semiparametric mixture cure models. *Computer Methods and
Programs in Biomedicine, 108*(3): 1255-1260.

-   For the implementation of the EM-algorithm, we build on previous
    work by:

Lai, X. & Yau, K. K. W. 2008. Long-term survivor model with bivariate
random effects: Applications to bone marrow transplant and carcinoma
study data. *Statistics in Medicine, 27*(27): 5692-5708.

and

Lai, X. & Yau, K. K. W. 2009. Multilevel Mixture Cure Models with Random
Effects. *Biometrical Journal, 51*(3): 456-466.

Lai, X. & Yau, K. K. W. 2010. Extending the long-term survivor mixture
model with random effects for clustered survival data. *Computational
Statistics and Data Analysis, 54*(9): 2103-2112.
