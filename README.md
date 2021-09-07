# Wilcoxmed package <img src="https://www.rdocumentation.org/badges/version/wilcoxmed">
The R package implementation of the 1-Sample Wilcoxon Signed Rank Hypothesis Test for medians - allows breaking of ties, computation of the Wilcoxon Test Statistic with normal approximations.

We refer to the Wilcoxon Sign Ranked hypothesis test for medians for the one-sample problem.

Given that <img src="https://render.githubusercontent.com/render/math?math=X_1, X_2, ... , X_n~F">.
Assume that  <img src="https://render.githubusercontent.com/render/math?math=F"> is unknown except that it is continuous and symmetric.
We wish to test  <img src="https://render.githubusercontent.com/render/math?math=H_0:m=m_0"> against <img src="https://render.githubusercontent.com/render/math?math=H_1:m>m_0">
at the significance level <img src="https://render.githubusercontent.com/render/math?math=\alpha">

Let <img src="https://render.githubusercontent.com/render/math?math=R_i=(signed \, rank \,  of X_i-m_0)">.
The test statistic is <img src="https://render.githubusercontent.com/render/math?math=\sum R_i">, 
or equivalently <img src="https://render.githubusercontent.com/render/math?math=W=\sum R_i/2%2Bn(n%2B1)/4">

We reject <img src="https://render.githubusercontent.com/render/math?math=H_0"> iff <img src="https://render.githubusercontent.com/render/math?math=p-value = P(W\geq w_{obs}|H_0)<\alpha">

For this problem, there could be ties in the data, however the algorithm is able to break the ties by assigning the average absolute rank to each tied value.
The exact distribution of <img src="https://render.githubusercontent.com/render/math?math=W"> is taken from Bickel and Doksum (1973).

The normal approximations to <img src="https://render.githubusercontent.com/render/math?math=P(W\leq w)"> can be found via the following theorem:

Theorem: If <img src="https://render.githubusercontent.com/render/math?math=X_1, X_2, ... , X_n~F">, <img src="https://render.githubusercontent.com/render/math?math=F"> is continuous,
then as <img src="https://render.githubusercontent.com/render/math?math=n\to\infty">, we have: 
<img src="https://render.githubusercontent.com/render/math?math=W\stackrel{approx}{\sim}N(n(n%2B1)/4 \, , \, n(n%2B1)(2n%2B1)/24)">

# Examples of code implementation

### Hypothesis Testing

Given some data:  <img src="https://render.githubusercontent.com/render/math?math=X\in\{3, 4, 7, 10, 4, 12, 1, 9, 2, 15\}"> 

If we want to test the hypotheses <img src="https://render.githubusercontent.com/render/math?math=H_0:m=5"> against
<img src="https://render.githubusercontent.com/render/math?math=H_1:m>5">
without using normal approximation:

vec = c(3, 4, 7, 10, 4, 12, 1, 9, 2, 15)

res = Wilcox.m.test(dat = vec, m_h0 = 5,
alternative = 'greater', normal_approx = F)

If we want to apply the normal approximation(Z-test), with the same hypotheses:

res = Wilcox.m.test(dat = vec, m_h0 = 5,
alternative = 'greater', normal_approx = T)

### To find the exact probabilities of the Wilcoxon sign distribution:

To find <img src="https://render.githubusercontent.com/render/math?math=P(W\leq3)"> for <img src="https://render.githubusercontent.com/render/math?math=n=5">:

W_stat(n=5, test_stat = 3, side = 'leq')


