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
devtools::install_github('HansTierens/MMixCure',dep=T)
```

    ## Downloading GitHub repo HansTierens/MMixCure@HEAD

    ## 
    ##          checking for file 'C:\Users\hans.tierens\AppData\Local\Temp\RtmpcxN1dI\remotesee384d7a3645\HansTierens-MMixCure-b7bc8d0/DESCRIPTION' ...  v  checking for file 'C:\Users\hans.tierens\AppData\Local\Temp\RtmpcxN1dI\remotesee384d7a3645\HansTierens-MMixCure-b7bc8d0/DESCRIPTION' (557ms)
    ##       -  preparing 'MMixCure':
    ##    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   v  checking DESCRIPTION meta-information
    ##       -  checking for LF line-endings in source and make files and shell scripts
    ##   -  checking for empty or unneeded directories
    ##       -  building 'MMixCure_0.1.0.tar.gz'
    ##      
    ## 

    ## Installing package into 'C:/Users/hans.tierens/Documents/R/win-library/4.0'
    ## (as 'lib' is unspecified)

``` r
library(MMixCure)
```

    ## Loading required package: arm

    ## Loading required package: MASS

    ## Loading required package: Matrix

    ## Loading required package: lme4

    ## 
    ## arm (Version 1.11-2, built: 2020-7-27)

    ## Working directory is C:/Users/hans.tierens/Desktop/MMixCure

    ## Loading required package: coxme

    ## Loading required package: survival

    ## Loading required package: bdsmatrix

    ## 
    ## Attaching package: 'bdsmatrix'

    ## The following object is masked from 'package:base':
    ## 
    ##     backsolve

    ## Loading required package: nloptr

    ## Loading required package: smcure

### Simulate some data

``` r
i=1
n=1000
ng=10
coefficients<-list(gamma0 = c(0,.40,-.30,-.30,.25),  beta0 = c(0,.45,.35,-.25,-.15))
revcovmat<-list(varcure=1,  varfrail=1,  corr=-.8)

data<-simsomdata(i,n,ng,coefficients,revcovmat)
knitr::kable(data[1:10,])
```

|          Y|  Delta|          Xs|          Xd|   Ds|   Dd|  grp|
|----------:|------:|-----------:|-----------:|----:|----:|----:|
|  0.2776793|      1|   0.8685774|   2.0591151|    1|    1|   10|
|  1.3678816|      1|  -1.0780345|  -0.0961369|    0|    1|    6|
|  0.4591765|      0|  -1.2223320|  -1.1254907|    1|    0|    4|
|  1.0618977|      0|  -0.7114447|   0.0502006|    0|    1|    4|
|  3.2618346|      0|  -1.4240317|  -0.1836796|    1|    0|   10|
|  1.2679387|      1|  -1.6693444|   0.2374040|    0|    0|    9|
|  0.1752351|      1|   1.3792361|  -0.5707300|    0|    0|    7|
|  0.6885131|      1|  -0.9196746|   0.6788652|    1|    1|    6|
|  0.2191299|      1|  -0.5044900|  -1.4653455|    1|    0|    9|
|  1.6844035|      1|  -1.1347318|  -1.1376878|    1|    1|    8|

``` r
knitr::kable(summary(data))
```

<table>
<colgroup>
<col style="width: 2%" />
<col style="width: 15%" />
<col style="width: 12%" />
<col style="width: 15%" />
<col style="width: 15%" />
<col style="width: 12%" />
<col style="width: 11%" />
<col style="width: 13%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: left;">Y</th>
<th style="text-align: left;">Delta</th>
<th style="text-align: left;">Xs</th>
<th style="text-align: left;">Xd</th>
<th style="text-align: left;">Ds</th>
<th style="text-align: left;">Dd</th>
<th style="text-align: left;">grp</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;"></td>
<td style="text-align: left;">Min. : 0.05813</td>
<td style="text-align: left;">Min. :0.000</td>
<td style="text-align: left;">Min. :-3.23639</td>
<td style="text-align: left;">Min. :-3.04536</td>
<td style="text-align: left;">Min. :0.000</td>
<td style="text-align: left;">Min. :0.0</td>
<td style="text-align: left;">Min. : 1.000</td>
</tr>
<tr class="even">
<td style="text-align: left;"></td>
<td style="text-align: left;">1st Qu.: 0.80589</td>
<td style="text-align: left;">1st Qu.:0.000</td>
<td style="text-align: left;">1st Qu.:-0.73986</td>
<td style="text-align: left;">1st Qu.:-0.60703</td>
<td style="text-align: left;">1st Qu.:0.000</td>
<td style="text-align: left;">1st Qu.:0.0</td>
<td style="text-align: left;">1st Qu.: 3.000</td>
</tr>
<tr class="odd">
<td style="text-align: left;"></td>
<td style="text-align: left;">Median : 1.38872</td>
<td style="text-align: left;">Median :0.000</td>
<td style="text-align: left;">Median :-0.02242</td>
<td style="text-align: left;">Median : 0.05012</td>
<td style="text-align: left;">Median :0.000</td>
<td style="text-align: left;">Median :0.5</td>
<td style="text-align: left;">Median : 6.000</td>
</tr>
<tr class="even">
<td style="text-align: left;"></td>
<td style="text-align: left;">Mean : 1.97544</td>
<td style="text-align: left;">Mean :0.486</td>
<td style="text-align: left;">Mean :-0.02956</td>
<td style="text-align: left;">Mean : 0.05570</td>
<td style="text-align: left;">Mean :0.488</td>
<td style="text-align: left;">Mean :0.5</td>
<td style="text-align: left;">Mean : 5.643</td>
</tr>
<tr class="odd">
<td style="text-align: left;"></td>
<td style="text-align: left;">3rd Qu.: 2.57744</td>
<td style="text-align: left;">3rd Qu.:1.000</td>
<td style="text-align: left;">3rd Qu.: 0.70423</td>
<td style="text-align: left;">3rd Qu.: 0.75665</td>
<td style="text-align: left;">3rd Qu.:1.000</td>
<td style="text-align: left;">3rd Qu.:1.0</td>
<td style="text-align: left;">3rd Qu.: 8.000</td>
</tr>
<tr class="even">
<td style="text-align: left;"></td>
<td style="text-align: left;">Max. :13.04814</td>
<td style="text-align: left;">Max. :1.000</td>
<td style="text-align: left;">Max. : 2.96174</td>
<td style="text-align: left;">Max. : 3.03903</td>
<td style="text-align: left;">Max. :1.000</td>
<td style="text-align: left;">Max. :1.0</td>
<td style="text-align: left;">Max. :10.000</td>
</tr>
</tbody>
</table>

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
    ##   n= 1000, number of events= 486 
    ## 
    ##        coef exp(coef) se(coef)      z Pr(>|z|)    
    ## Xs  0.40004   1.49188  0.04597  8.702  < 2e-16 ***
    ## Xd  0.00425   1.00426  0.04416  0.096    0.923    
    ## Ds -0.40327   0.66813  0.09184 -4.391 1.13e-05 ***
    ## Dd  0.05369   1.05516  0.09080  0.591    0.554    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##    exp(coef) exp(-coef) lower .95 upper .95
    ## Xs    1.4919     0.6703    1.3633    1.6325
    ## Xd    1.0043     0.9958    0.9210    1.0950
    ## Ds    0.6681     1.4967    0.5581    0.7999
    ## Dd    1.0552     0.9477    0.8831    1.2607
    ## 
    ## Concordance= 0.64  (se = 0.013 )
    ## Likelihood ratio test= 95.33  on 4 df,   p=<2e-16
    ## Wald test            = 95.31  on 4 df,   p=<2e-16
    ## Score (logrank) test = 95.48  on 4 df,   p=<2e-16

``` r
frailtyfit<-coxph(Surv(Y,Delta)~1+Xs+Xd+Ds+Dd+frailty.gaussian(grp),data=data,ties='breslow')
summary(frailtyfit)
```

    ## Call:
    ## coxph(formula = Surv(Y, Delta) ~ 1 + Xs + Xd + Ds + Dd + frailty.gaussian(grp), 
    ##     data = data, ties = "breslow")
    ## 
    ##   n= 1000, number of events= 486 
    ## 
    ##                       coef    se(coef) se2     Chisq DF   p      
    ## Xs                     0.4294 0.04705  0.04700 83.29 1.00 7.1e-20
    ## Xd                     0.0143 0.04443  0.04438  0.10 1.00 7.5e-01
    ## Ds                    -0.3856 0.09283  0.09270 17.25 1.00 3.3e-05
    ## Dd                     0.0692 0.09121  0.09116  0.58 1.00 4.5e-01
    ## frailty.gaussian(grp)                          63.98 8.66 1.6e-10
    ## 
    ##    exp(coef) exp(-coef) lower .95 upper .95
    ## Xs    1.5363     0.6509    1.4009    1.6847
    ## Xd    1.0144     0.9858    0.9298    1.1067
    ## Ds    0.6801     1.4705    0.5669    0.8158
    ## Dd    1.0716     0.9331    0.8962    1.2814
    ## 
    ## Iterations: 5 outer, 18 Newton-Raphson
    ##      Variance of random effect= 0.1736209 
    ## Degrees of freedom for terms= 1.0 1.0 1.0 1.0 8.7 
    ## Concordance= 0.669  (se = 0.669 )
    ## Likelihood ratio test= 173.8  on 12.66 df,   p=<2e-16

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
    ##               Estimate  Std.Error   Z value     Pr(>|Z|)
    ## (Intercept)  0.3978187 0.12933041  3.075987 2.098067e-03
    ## Xs           0.4441501 0.07353542  6.039947 1.541649e-09
    ## Xd          -0.1800797 0.06819397 -2.640698 8.273534e-03
    ## Ds          -0.4646626 0.16775771 -2.769843 5.608325e-03
    ## Dd           0.1566690 0.13840869  1.131931 2.576636e-01
    ## 
    ## 
    ## Failure time distribution model:
    ##        Estimate  Std.Error     Z value     Pr(>|Z|)
    ## Xs  0.349056245 0.05372259  6.49738315 8.172907e-11
    ## Xd  0.302361388 0.05617566  5.38242705 7.348816e-08
    ## Ds -0.314271919 0.09776829 -3.21445660 1.306917e-03
    ## Dd -0.009106158 0.09229785 -0.09866056 9.214078e-01

``` r
printsmcure(spmcfit)
```

    ## Call:
    ## smcure(formula = Surv(Y, Delta) ~ 1 + Xs + Xd + Ds + Dd, cureform = ~1 + 
    ##     Xs + Xd + Ds + Dd, data = data, model = "ph", link = "logit", 
    ##     Var = TRUE)
    ## 
    ## Cure probability model:
    ##               Estimate  Std.Error   Z value     Pr(>|Z|)
    ## (Intercept)  0.3978187 0.12933041  3.075987 2.098067e-03
    ## Xs           0.4441501 0.07353542  6.039947 1.541649e-09
    ## Xd          -0.1800797 0.06819397 -2.640698 8.273534e-03
    ## Ds          -0.4646626 0.16775771 -2.769843 5.608325e-03
    ## Dd           0.1566690 0.13840869  1.131931 2.576636e-01
    ## 
    ## 
    ## Failure time distribution model:
    ##        Estimate  Std.Error     Z value     Pr(>|Z|)
    ## Xs  0.349056245 0.05372259  6.49738315 8.172907e-11
    ## Xd  0.302361388 0.05617566  5.38242705 7.348816e-08
    ## Ds -0.314271919 0.09776829 -3.21445660 1.306917e-03
    ## Dd -0.009106158 0.09229785 -0.09866056 9.214078e-01

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
    ## (Intercept)  0.09973   0.32672  0.30526  0.76017
    ## Xs           0.47928   0.07143  6.70997  0.00000
    ## Xd          -0.17574   0.07210 -2.43732  0.01480
    ## Ds          -0.45999   0.14111 -3.25986  0.00111
    ## Dd           0.14420   0.14107  1.02223  0.30667
    ## 
    ## 
    ## Failure time distribution model:
    ##    Estimate Std.Error  Z value Pr(>|Z|)
    ## Xs  0.40283   0.04849  8.30730  0.00000
    ## Xd  0.34978   0.04846  7.21860  0.00000
    ## Ds -0.38821   0.09493 -4.08965  0.00004
    ## Dd -0.11702   0.09328 -1.25444  0.20968
    ## 
    ## Random Effects:
    ##                  Variance Correlation
    ## Cure/Incidence      0.917            
    ## Survival/Latency    0.960      -0.838

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
