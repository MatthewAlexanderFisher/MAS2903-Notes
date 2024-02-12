# Chapter 2: Bayes Theorem for Distributions

## Section 2.1: Introduction

Suppose we have data $\underline{x}$ which we model using the
probability (density) function $f(\underline{x}|\theta)$, which depends
on a single parameter $\theta$. Once we have observed the data,
$f(\underline{x}|\theta)$ is the *likelihood function* for $\theta$ and
is a function of $\theta$ (for fixed $\underline{x}$) rather than
of $\underline{x}$ (for fixed $\theta$).

Also, suppose we have prior beliefs about likely values of $\theta$
expressed by a probability (density) function $\pi(\theta)$. We can
combine both pieces of information using the following version of Bayes
Theorem. The resulting distribution for $\theta$ is called the posterior
distribution for $\theta$ as it expresses our beliefs about $\theta$
*after* seeing the data. It summarises all our current knowledge about
the parameter $\theta$.

### Bayes Theorem 

The posterior probability (density) function for $\theta$ is

$$\pi(\theta|\underline{x})=\frac{\pi(\theta)\,f(\underline{x}|\theta)}{f(\underline{x})}$$

where 

$$
f(\underline{x})=
\begin{cases}
\int_\Theta\pi(\theta)\,f(\underline{x}|\theta)\,d\theta & \text{if
$\theta$ is continuous}, \\ \\
\sum_\Theta\pi(\theta)\,f(\underline{x}|\theta) & \text{if $\theta$ is discrete}.
\end{cases}$$

 Notice that, as $f(\underline{x})$ is not a function of
$\theta$, Bayes Theorem can be rewritten as 

$$
\begin{aligned}
\pi(\theta|\underline{x})&\propto \pi(\theta)\times f(\underline{x}|\theta) \\
i.e. \, \text{posterior}&\propto\text{prior}\times\text{likelihood}.
\end{aligned}$$

Thus, to obtain the posterior distribution, we need:

-   data, from which we can form the *likelihood*
    $f(\underline{x}|\theta)$, and

-   a suitable distribution, $\pi(\theta)$, that represents our *prior
    beliefs* about $\theta$.

You should now be comfortable with how to obtain the likelihood. But how
do we specify a prior (point 2)? In Chapter 3 we will consider different
approaches to specifying prior distributions.

For now, we will assume someone else has done this for us; the main aim
of this chapter is simply to operate Bayes Theorem for distributions to
obtain the posterior distribution for $\theta$.

Before we do this, it will be worth re--familiarising ourselves with
some continuous probability distributions you have met before, and which
we will use extensively in this course: the uniform, beta and gamma
distributions (indeed, I will assume that you are more than familiar
with some other 'standard' distributions we will use -- e.g. the
exponential, Normal, Poisson, and binomial, and so will *not* review
these here).

##### Definition 2.1: Continuous Uniform distribution

 \
The random variable $Y$ follows a Uniform distribution, denoted
$U(a,b)$, if it has probability density function

$$
f(y|a,b) = \frac{1}{b-a}, \quad \quad a\leq y \leq b.$$

 This form of
probability density function ensures that all values in the range
$[a,b]$ are *equally likely*, hence the name "uniform". This
distribution is sometimes called the *rectangular distribution* because
of its shape.

You should remember from MAS1608 that

$$
\text{E}[Y] = \frac{a+b}{2} \quad \mathrm{and} \quad \text{Var}[Y] = \frac{(b-a)^{2}}{12}.$$

In the space below, sketch the probability density functions for
$U(0,1)$ and $U(10,50)$.

::: gapbox
:::

##### Definition 2.2: Beta distribution

 \
The random variable $Y$ follows a beta $\mathrm{Beta}(a,b)$ distribution
($a>0$,$b>0$) if it has probability density function

```{math}
:label: eq:betapdf
f(y|a,b)=\frac{y^{a-1}(1-y)^{b-1}}{\mathrm{B}(a,b)},\quad\quad 0<y<1.
```

 The constant term
$\mathrm{B}(a,b)$, also known as the *beta function*, ensures that the
density integrates to one. Therefore 

```{math}
:label: eq:betafn
\mathrm{B}(a,b)=\int_0^1
y^{a-1}(1-y)^{b-1}\,dy.
```

 It can be shown that the beta
function can be expressed in terms of another function, called the
*gamma function* $\Gamma(\cdot)$, as

$$
\mathrm{B}(a,b) =\frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)},$$

 where

```{math}
:label: eq:gammafn
\Gamma(a)=\int_0^\infty x^{a-1}e^{-x}\,dx.
```

 Tables are available for
both $\mathrm{B}(a,b)$ and $\Gamma(a)$. However, these functions are
very simple to evaluate when $a$ and $b$ are integers since the gamma
function is a generalisation of the factorial function. In particular,
when $a$ and $b$ are integers, we have

$$
\Gamma(a)=(a-1)!\quad\text{and}\quad 
\mathrm{B}(a,b)=\frac{(a-1)!(b-1)!}{(a+b-1)!}.$$

 For example,

$$
\mathrm{B}(2,3)=\frac{1!\times2!}{4!}=\frac{1}{12}.$$

 It can be shown,
using the identity $\Gamma(a)=(a-1)\Gamma(a-1)$, that 

$$
\begin{aligned}
\text{E}[Y]&=\frac{a}{a+b},\quad\quad\text{and}\quad\quad 
\text{Var}[Y]=\frac{ab}{(a+b)^2(a+b+1)}. 
\end{aligned}$$

 Also,

$$
\textnormal{Mode}(Y)=\frac{a-1}{a+b-2}, \quad\text{if $a>1$ and $b>1$}.$$

##### Definition 2.3: Gamma distribution

 \
The random variable $Y$ follows a Gamma $\mathrm{Gamma}(a,b)$
distribution ($a>0$, $b>0$) if it has probability density function

$$
f(y|a,b) =\frac{b^ay^{a-1}e^{-by}}{\Gamma(a)},\quad\quad y>0,$$

 where
$\Gamma(a)$ is the gamma function defined in . It can be shown that

$$
\begin{aligned}
\text{E}[Y]&=\frac{a}{b}\quad\quad\text{and}\quad\quad 
\text{Var}[Y]=\frac{a}{b^2}. 
\end{aligned}$$

 Also,

$$
\textnormal{Mode}(Y)=\frac{a-1}{b}, \quad\text{if $a\geq 1$}.$$

We can use `R` to visualise the beta and gamma distributions for various
values of $(a,b)$ (and indeed any other standard probability
distribution you have met so far). For example, we know that the beta
distribution is valid for all values in the range $(0,1)$. In `R`, we
can set this up by typing:

``` Python
x = seq(0,1, 0.01)
```

which specifies $x$ to take all values in the range 0 to 1, in steps of
0.01. The following code then calculates the density of
$\mathrm{Beta}(2,5)$, as given by Equation with $a=2$ and $b=5$:

``` Python
y = dbeta(x, 2,5)
```

Plotting $y$ against $x$ and joining with lines gives the
$\mathrm{Beta}(2,5)$ density shown in Figure 2.1 (top left); in `R` this
is achieved by typing:

``` Python
plot(x, y, type='l')
```

Also shown in Figure 2.1 are densities for the $\mathrm{Beta}(0.5,0.5)$
(top right), $\mathrm{Beta}(77,5)$ (bottom left) and
$\mathrm{Beta}(10,10)$ (bottom right) distributions. Notice that
different combinations of $(a,b)$ give rise to different shapes of
distribution between the limits of 0 and 1 -- symmetric, positively
skewed and negatively skewed: careful choices of $a$ and $b$ could thus
be used to express our prior beliefs about probabilities/proportions we
think might be more or less likely to occur. When $a=b$ we have a
distribution which is symmetric about 0.5. Similar plots can be
constructed for any standard distribution of interest using, for
example, `dgamma` or `dnorm` instead of `dbeta` for the gamma or Normal
distributions, respectively; Figure 2.2 shows densities for various
gamma distributions.

![Plots of $\mathrm{Beta}(a,b)$ densities for various values of
$(a,b)$.](images/figures/betaplot1.svg)

![Plots of $\mathrm{Gamma}(a,2)$ densities, for various values of
$a$.](images/gammaplot1.svg)

## Section 2.2: Bayes Theorem for distributions in action

We will now see Bayes Theorem for distributions in operation. Remember
-- for now, we will assume that someone else has provided the prior
distribution for $\theta$. In Chapter 3 we will consider how this might
be done.

##### Example 2.1

 \
Consider an experiment with a possibly biased coin. Let

$\theta=\text{Pr}(\mathrm{Head})$. Suppose that, before conducting the
experiment, we believe that all values of $\theta$ are equally likely:
this gives a prior distribution $\theta\sim U(0,1)$, and so

```{math}
:label: eq:p1
\pi(\theta)=1,\quad\quad0<\theta<1.
```

 Note that with this prior
distribution $\text{E}[\theta]=0.5$. We now toss the coin 5 times and
observe 1 head. Determine the posterior distribution for $\theta$ given
this data.

```{dropdown} Solution

The data is an observation on the random variable $X|\theta\sim         \mathrm{Bin}(5,\theta)$. This gives a likelihood function

<div style="text-align: center;">
    <a id="eq:p2"></a>
    <div class="math notranslate nohighlight">
        $$
f(x=1|\theta)=5\theta(1-\theta)^4
        $$
    </div>
</div>

 which favours values of
$\theta$ near its maximum $\theta=0.2$. Therefore, we have a conflict of
opinions: the prior distribution suggests that $\theta$ is probably
around 0.5 and the data suggest that it is around 0.2. We can use Bayes
Theorem to combine these two sources of information in a coherent way.
First 

$$
\begin{aligned}
        f(x=1) &=\int_\Theta\pi(\theta)f(x=1\, |\,\theta)\,d\theta = \int_0^1 1\times 5\theta(1-\theta)^4\,d\theta \\
        &=\int_0^1 \theta\times 5(1-\theta)^4\,d\theta =\left[-(1-\theta)^5\,\theta\right]^1_0
        +\int_0^1 (1-\theta)^5\,d\theta \\
        &=0 + \left[-\frac{(1-\theta)^6}{6}\right]^1_0 =\frac{1}{6}. 
        
\end{aligned}$$

 Therefore, the posterior density is (for $0<\theta<1$):

$$
\begin{aligned}
        \pi(\theta|x=1)&=\frac{\pi(\theta)f(x=1|\theta)}{f(x=1)} =\frac{5\theta(1-\theta)^4}{1/6}\\
        &=30\,\theta(1-\theta)^4 =\frac{\theta(1-\theta)^4}{\mathrm{B}(2,5)},\quad\quad 0<\theta<1,
        
\end{aligned}$$

 and so the posterior distribution is
$\theta|x=1\sim \mathrm{Beta}(2,5)$ -- see Definition . This
distribution has its mode at $\theta=0.2$, and mean at
$\text{E}[\theta|x=1]=2/7=0.286$.

```

The main difficulty in calculating the posterior distribution was in
obtaining the $f(x)$ term. However, in many cases we can recognise the
posterior distribution without the need to calculate this constant term
(constant with respect to $\theta$). In this example, we can calculate
the posterior distribution as 

$$
\begin{aligned}
\pi(\theta|\underline{x})&\propto\pi(\theta)f(x=1|\theta) \\
&\propto 1\times 5\theta(1-\theta)^4,\quad\quad 0<\theta<1  \\
&=k\theta(1-\theta)^4,\quad\quad 0<\theta<1.
\end{aligned}$$

 As $\theta$ is a continuous quantity, what we would like
to know is what continuous distribution defined on $(0,1)$ has a
probability density function which takes the form
$k\theta^{g-1}(1-\theta)^{h-1}$. The answer is the $\mathrm{Beta}(g,h)$
distribution. Therefore, choosing $g$ and $h$ appropriately, we can see
that the posterior distribution is $\theta|x=1\sim \mathrm{Beta}(2,5)$.

**Summary:**

It is possible that we have a biased coin. If we suppose that all values
of $\theta=\text{Pr(Head)}$ are equally likely and then observe 1 head
out of 5, then the most likely value of $\theta$ is 0.2 --- the same as
the most likely value from the data alone (not surprising!). However, on
average, we would expect $\theta$ to be around 0.286. Uncertainty about
$\theta$ has changed from a (prior) standard deviation of 0.289 to a
(posterior) standard deviation of 0.160. The changes in our beliefs
about $\theta$ are more fully described by the prior and posterior
distributions shown in Figure .

![Prior (dashed) and posterior (solid) densities for
$\theta=\text{Pr(Head)}$](images/figures/priorplot1.svg)

##### Example 2.2

 \
Consider an experiment to determine how good a music expert is at
distinguishing between pages from Haydn and Mozart scores. Let
$\theta=\text{Pr}(\text{correct choice})$. Suppose that, before
conducting the experiment, we have been told that the expert is very
competent. In fact, it is suggested that we should have a prior
distribution which has a mode around $\theta=0.95$ and for which
$\textnormal{Pr}(\theta<0.8)$ is very small. We choose
$\theta\sim \mathrm{Beta}(77,5)$, with probability density function

```{math}
:label: eq:p3
\pi(\theta)=128107980\,\theta^{76}(1-\theta)^4,\quad\quad 0<\theta<1.
```

 A graph of this prior density is
given in Figure  .

![Prior density for the music expert's skill.](images/priorplot2.svg)

In the experiment, the music expert makes the correct choice 9 out of 10
times. Determine the posterior distribution for $\theta$ given this
information.

```{dropdown} Solution

We have an observation on the random variable $X|\theta\sim         \mathrm{Bin}(10,\theta)$. This gives a likelihood function of

<div style="text-align: center;">
    <a id="eq:p4"></a>
    <div class="math notranslate nohighlight">
        $$
f(x=9|\theta)=10\,\theta^9(1-\theta)
        $$
    </div>
</div>

 which favours values of
$\theta$ near its maximum $\theta=0.9$. We combine these two sources of
information using Bayes Theorem. The posterior density function is

<div style="text-align: center;">
    <a id="eq:p5"></a>
    <div class="math notranslate nohighlight">
        $$
\begin{aligned}
        \pi(\theta|x=9)&\propto\pi(\theta)f(x=9|\theta) \notag \\
        &\propto 128107980\,\theta^{76}(1-\theta)^4\times 10\,\theta^9(1-\theta),
        \quad\quad 0<\theta<1 \notag \\
        &=k\theta^{85}(1-\theta)^5,\quad\quad 0<\theta<1.
\end{aligned}
        $$
    </div>
</div>

 We can recognise this density function as one from the
Beta family. So, the posterior distribution is $\theta|x=9\sim         \mathrm{Beta}(86,6)$.

```

**Summary:**

The changes in our beliefs about $\theta$ are described by the prior and
posterior distributions shown in Figure  and summarised in Table .

![Prior (dashed) and posterior (solid) densities for the music expert's
skill for the $0.8 < \theta < 1$.](images/priorposterior1.svg)

```{table} Changes in beliefs about $\theta$.
:widths: auto
:align: center
:name: tab:betaplot4

|                           | [Prior](#eq:p3) | [Likelihood](#eq:p4) | [Posterior](#eq:p5) |
|---------------------------|---------------------|--------------------------|-------------------------|
| $\textnormal{Mode}(\theta)$ | 0.950              | 0.900                    | 0.944                   |
| $\textnormal{E}[\theta]$              | 0.939               | --                       | 0.935                   |
| $\textnormal{SD}(\theta)$ | 0.0263              | --                       | 0.0256                  |
```


Notice that, having observed only a 90% success rate in the experiment,
the posterior mode and mean are smaller than their prior values. Also,
the experiment has largely confirmed our ideas about $\theta$, with the
uncertainty about $\theta$ being only very slightly reduced.

##### Example 2.3

 \
Max, a video game pirate, is trying to identify the proportion of
potential customers $\theta$ who might be interested in buying *Call of
Duty: Modern Warfare II* next month. Based on the proportion of
customers who have bought similarly violent games from him in the past,
he assumes that $\theta \sim \mathrm{Beta}(2.5,12)$; a plot of this
prior density is shown in Figure .

![Max's prior density.](images/priorplot3.svg)

Max asks five potential customers if they would buy *Call of Duty:
Modern Warfare II* from him, and four say they would. Using this
information, what is Max's posterior distribution for $\theta$?

```{dropdown} Solution

We have been told that the prior for $\theta$ is a
$\mathrm{Beta}(2.5,12)$ distribution -- this has density given by

<div style="text-align: center;">
    <a id="eq:max1"></a>
    <div class="math notranslate nohighlight">
        $$
\frac{1}{\mathrm{B}(2.5,12)}\theta^{2.5-1}(1-\theta)^{12-1} = 435.1867\theta^{1.5}(1-\theta)^{11}.
        $$
    </div>
</div>

 We have an observation on the
random variable $X|\theta \sim \mathrm{Bin}(5,\theta)$. This gives a
likelihood function of

<div style="text-align: center;">
    <a id="eq:max2"></a>
    <div class="math notranslate nohighlight">
        $$
f(x=4|\theta)  = \binom{5}{4} \theta^{4}(1-\theta)^{1} = 5 \theta^{4}(1-\theta),
        $$
    </div>
</div>

 which favours values of
$\theta$ near its maximum 0.8. We combine our prior information () with
the data () -- to obtain our posterior distribution -- using Bayes
Theorem. The posterior density function is 

$$
\begin{aligned}
    \pi(\theta|x=4) &\propto& \pi(\theta)f(x=4|\theta)\\
                    &\propto& 435.1867 \theta^{1.5}(1-\theta)^{11} \times 5 \theta^{4}(1-\theta), \qquad 0<\theta<1,  \\
    
\end{aligned}$$

 giving

<div style="text-align: center;">
    <a id="eq:max3"></a>
    <div class="math notranslate nohighlight">
        $$
\pi(\theta|x=4) = k \theta^{5.5}(1-\theta)^{12}, \qquad 0<\theta<1.
        $$
    </div>
</div>

 You should recognise this
density function as one from the beta family. In fact, we have a
$\mathrm{Beta}(6.5, 13)$, i.e. $\theta|x=4 \sim \mathrm{Beta}(6.5,13)$.

```

**Summary:**

The changes in our beliefs about $\theta$ are described by the prior and
posterior distributions shown in Figure 2.7 and summarised in Table 2.2.

![Prior (dashed) and posterior (solid) densities for Max's
problem.](images/priorposterior2.svg)

```{table} Changes in beliefs about $\theta$.
:widths: auto
:align: center
:name: tab:betaposterior

|                           | [Prior](#eq:max1) | [Likelihood](#eq:max2) | [Posterior](#eq:max3) |
|---------------------------|---------------------|--------------------------|-------------------------|
| $\textnormal{Mode}(\theta)$ | 0.12               | 0.8                      | 0.314                   |
| $\textnormal{E}[\theta]$              | 0.172               | --                       | 0.333                   |
| $\textnormal{SD}(\theta)$ | 0.096               | --                       | 0.104                   |
```


Notice how the posterior has been "pulled" from the prior towards the
observed value: the mode has moved up from 0.12 to 0.314, and the mean
has moved up from 0.172 to 0.333. Having just one observation in the
likelihood, we see that there is hardly any change in the standard
deviation from prior to posterior: we would expect to see a decrease in
standard deviation with the addition of more data values.

##### Example 2.4

 \
Table  shows some data on the times between serious earthquakes. An
earthquake is included if its magnitude is at least 7.5 on the Richter
scale or if over 1000 people were killed. Recording starts on
16 December 1902 (4500 killed in Turkistan). The table includes data on
21 earthquakes, that is, 20 "waiting times" between earthquakes.

```{table} Time intervals between major earthquakes (in days).
:widths: auto
:align: center
:name: tab:earth

| 840 | 157 | 145 | 44  | 33  | 121 | 150 | 280 | 434 | 736 |
|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| 584 | 887 | 263 | 1901| 695 | 294 | 562 | 721 | 76  | 710 |
```


It is believed that earthquakes happen in a random haphazard kind of way
and that times between earthquakes can be described by an exponential
distribution. Data over a much longer period suggest that this
exponential assumption is plausible. Therefore, we will assume that
these data are a random sample from an exponential distribution with
rate $\theta$ (and mean $1/\theta$). The parameter $\theta$ describes
the rate at which earthquakes occur.

An expert on earthquakes has prior beliefs about the rate of
earthquakes, $\theta$, described by a $\mathrm{Gamma}(10,4000)$
distribution, which has density density

```{math}
:label: eq:earthprior
\pi(\theta)=\frac{4000^{10}\,\theta^9e^{-4000\theta}}{\Gamma(10)}, 
\quad\theta>0,
```

 and mean $\text{E}[\theta]=0.0025$. A plot of this
prior distribution can be found in Figure . As you might expect, the
expert believes that, realistically, only very small values of $\theta$
are likely, though larger values are not ruled out! Determine the
posterior distribution for $\theta$.

![Prior density for the earthquake rate $\theta$](images/priorplot4.svg)

```{dropdown} Solution

The data are observations on $X_i|\theta\sim Exp(\theta)$,
$i=1,2,\ldots,20$ (independent). Therefore, the likelihood function
for $\theta$ is 

<div style="text-align: center;">
    <a id="eq:earthlik"></a>
    <div class="math notranslate nohighlight">
        $$
\begin{aligned}
    f(\underline{x}|\theta)&=\prod_{i=1}^{20} \theta e^{-\theta x_i}, 
    \quad\quad\theta>0 \notag \\
    &=\theta^{20}\exp\left(-\theta\sum_{i=1}^{20} x_i\right),
    \quad\quad\theta>0 \notag \\
    &=\theta^{20} e^{-9633\theta},
    \quad\quad\theta>0. 
\end{aligned}
        $$
    </div>
</div>

 We now apply Bayes Theorem to combine the expert opinion
with the observed data. The posterior density function is

<div style="text-align: center;">
    <a id="eq:earthpost"></a>
    <div class="math notranslate nohighlight">
        $$
\begin{aligned}
    \pi(\theta|\underline{x})&\propto\pi(\theta)f(\underline{x}|\theta) \notag \\
    &\propto\frac{4000^{10}\,\theta^{9}e^{-4000\theta}}{\Gamma(10)}\times
    \theta^{20} e^{-9633\theta},\quad\quad\theta>0 \notag  \\
    &=k\,\theta^{30-1}e^{-13633\theta},\quad\quad\theta>0.
\end{aligned}
        $$
    </div>
</div>

 The only continuous distribution which takes the form
$k\theta^{g-1}e^{-h\theta}$, $\theta>0$ is the $\mathrm{Gamma}(g,h)$
distribution. Therefore, the posterior distribution must be
$\theta|\underline{x}\sim     \mathrm{Gamma}(30,13633)$.

```

**Summary:**

The data have updated our beliefs about $\theta$ from a
$\mathrm{Gamma}(10,4000)$ distribution to a\
$\mathrm{Gamma}(30,13633)$ distribution. Plots of these distributions
are given in Figure , and Table  gives a summary of the main changes
induced by incorporating the data. Notice that, as the mode of the
likelihood function is close to that of the prior distribution, the
information in the data is consistent with that in the prior
distribution. This results in a reduction in variability from the prior
to the posterior distributions. The similarity between the prior beliefs
and the data has reduced the uncertainty we have about the likely
earthquake rate $\theta$.

![Prior (dashed) and posterior (solid) densities for the earthquake rate
$\theta$.](images/priorposterior3.svg)

```{table} Changes in beliefs about $\theta$.
:widths: auto
:align: center
:name: tab:earthsum

|                           | [Prior](#eq:earthprior) | [Likelihood](#eq:earthlik) | [Posterior](#eq:earthpost) |
|---------------------------|---------------------------|------------------------------|------------------------------|
| $\textnormal{Mode}(\theta)$ | 0.00225                   | 0.00208                      | 0.00213                      |
| $\textnormal{E}[\theta]$              | 0.00250                   | --                           | 0.00220                      |
| $\textnormal{SD}(\theta)$ | 0.00079                   | --                           | 0.00040                      |
```


##### Example 2.5

 \
We now consider the general case of the problem discussed in Example .
Suppose $X_i|\theta\sim \mathrm{Exp}(\theta)$, $i=1,2,\ldots,n$
(independent) and our prior beliefs about $\theta$ are summarised by a
$\mathrm{Gamma}(g,h)$ distribution (with $g$ and $h$ known), with
density 

```{math}
:label: eq:p6
\pi(\theta)=\frac{h^g\,\theta^{g-1}e^{-h\theta}}{\Gamma(g)}, 
\quad\theta>0.
```

 Determine the posterior distribution for $\theta$.

```{dropdown} Solution

The likelihood function for $\theta$ is 

<div style="text-align: center;">
    <a id="eq:p7"></a>
    <div class="math notranslate nohighlight">
        $$
\begin{aligned}
    f(\underline{x}|\theta)&=\prod_{i=1}^n \theta e^{-\theta x_i}, 
    \quad\quad\theta>0 \notag \\
    &= \theta^n e^{-n\bar x\theta}, \quad\quad\theta>0. 
\end{aligned}
        $$
    </div>
</div>

 We now apply Bayes Theorem. The posterior density
function is 

$$
\begin{aligned}
    \pi(\theta|\underline{x})&\propto\pi(\theta)f(\underline{x}|\theta) \notag \\
    &\propto\frac{h^g\,\theta^{g-1}e^{-h\theta}}{\Gamma(g)}\times
    \theta^n e^{-n\bar x\theta},\quad\quad\theta>0. \notag
    
\end{aligned}$$

 Therefore,

<div style="text-align: center;">
    <a id="eq:p8"></a>
    <div class="math notranslate nohighlight">
        $$
\pi(\theta|\underline{x}) = k\theta^{g+n-1}e^{-(h+n\bar x)\theta},\quad\quad\theta>0
        $$
    </div>
</div>

where $k$ is a constant that does not depend on $\theta$. Therefore, the
posterior distribution takes the form $k\theta^{g-1}e^{-h\theta}$,
$\theta>0$ and so must be a gamma distribution. Thus we have
$\theta|\underline{x}\sim \mathrm{Gamma}(g+n,h+n\bar x)$.

```

**Summary:**

If we have a random sample from an $\mathrm{Exp}(\theta)$ distribution
and our prior beliefs about $\theta$ follow a $\mathrm{Gamma}(g,h)$
distribution then, after incorporating the data, our (posterior) beliefs
about $\theta$ follow a $\mathrm{Gamma}(g+n,h+n\bar x)$ distribution.
The changes in our beliefs about $\theta$ are summarised in Table ,
taking $g\geq 1$.

```{table} Changes in beliefs about $\theta$.
:widths: auto
:align: center
:name: tab:gam

|                           | [Prior](#eq:p6) | [Likelihood](#eq:p7) | [Posterior](#eq:p8) |
|---------------------------|--------------------|------------------------|-----------------------|
| $\textnormal{Mode}(\theta)$ | $(g-1)/h$           | $1/\bar x$             | $(g+n-1)/(h+n\bar x)$ |
| $\textnormal{E}[\theta]$              | $g/h$               | --                     | $(g+n)/(h+n\bar x)$   |
| $\textnormal{SD}(\theta)$ | $\sqrt{g}/h$        | --                     | $\sqrt{g+n}/(h+n\bar x)$ |
```


Notice that the posterior mean is greater than the prior mean if and
only if the likelihood mode is greater than the prior mean, that is,

$$
\text{E}[\theta|\underline{x}]>\text{E}[\theta]\quad\iff\quad 
\textnormal{Mode}(f(\underline{x}|\theta))>\text{E}[\theta].$$

 The
standard deviation of the posterior distribution is smaller than that of
the prior distribution if and only if the sample mean is large enough,
that is

$$
\textnormal{SD}(\theta|\underline{x})<\textnormal{SD}(\theta)\quad\iff\quad \bar x > k.$$

##### Example 2.6

 \
Software engineers at **Twitter** are interested in the click-through
rate $\theta$ of their advertisement recommender.

The **click-through rate** is the proportion of clicks out of the total
number of impressions of an online advertisement.

The engineers assess the click-through rate by testing their recommender
on users of Twitter. They obtain a total of $n$ impressions.

-   Using a **Beta prior**, quantify **your prior beliefs** about
    $\theta$.

```{dropdown} Solution

By plotting different Beta distributions, or by specifying
$\text{E}[\theta]$ and $\text{Var}[\theta]$ and solving the system of
equations to specify $a$ and $b$, I arrived at the following prior:

$$
\theta \sim \mathrm{Beta}(5, 100).$$

 This places the vast majority of
the probability mass for values $\theta < 0.2$. I chose this because I
expect that the click-through rate is very small, but I also have some
uncertainty about $\theta$.

```

-   Specify a **probability model** for the number of successful clicks.

```{dropdown} Solution

In each given impression, a user either clicks on the advert or does
not. A success is if the user does click on the advert. Therefore,
assuming independence between each impression (which is plausible), the
appropriate model is 

$$
X|\theta \sim \mathrm{Binomial}(n, \theta),$$

where $n$ is the total number of impressions.

```

-   Using your probability model, and a general
    $\theta \sim \mathrm{Beta}(a,b)$ prior, derive the posterior
    distribution.

```{dropdown} Solution

The prior is $\theta\sim \mathrm{Beta}(a,b)$ and so the prior pdf is

$$
\pi(\theta) = \frac{\theta^{a-1}(1-\theta)^{b-1}}{\mathrm{B}(a,b)}.$$

Our probability model was $X|\theta \sim \mathrm{Binomial}(n, \theta)$.
Therefore, the posterior is 

$$
\begin{aligned}
        \pi(\theta|x) &\propto \pi(\theta) f(x|\theta) \\
        &\propto \frac{\theta^{a-1}(1-\theta)^{b-1}}{\mathrm{B}(a,b)} \times \binom{n}{x} \theta^{x}(1-\theta)^{n-x} \\
        &\propto \theta^{a + x - 1} (1-\theta)^{b + n - x - 1}.
    
\end{aligned}$$

 We recognise this as the same form as the Beta pdf and
we conclude that 

$$
\theta | x \sim \mathrm{Beta}(a + x, b + n - x).$$

```

-   Using your prior and the following data, derive your **posterior
    distribution**. Out $n=10,000$ impressions, there were 47 successful
    clicks.

```{dropdown} Solution

My prior was $\theta \sim \mathrm{Beta}(5, 100)$. Using the update
formula from part (c), with $n = 10000$ and $x = 47$, my revised belief
about $\theta$ is 

$$
\theta | x \sim \mathrm{Beta}(52, 10053).$$

 Note
that $\text{E}[\theta|x] = \frac{52}{10053} = 0.00517$ (3 s.f.) and
$\text{Var}[\theta|x] = 5.07\times 10^{-7}$ (3 s.f.). Due to the amount
of data, the posterior is dominated by the likelihood function since the
variance is so small.

```

**Summary:**

The prior, likelihood function and posterior can be seen in Figure .
Here, I used my prior $\theta \sim \mathrm{Beta}(5,100)$. Our changes in
belief about the click-through rate $\theta$ are summarised in Table ,
after incorporating the data.

Note how our posterior beliefs are dominated by the likelihood function.
This is because we have so much data ($n = 10000$).

![Prior, likelihood function and posterior of the click-through rate
$\theta$.](images/clickthroughrate_priorposterior.svg)

```{table} Changes in beliefs about $\theta$.
:widths: auto
:align: center
:name: tab:clickthrough

|                           | [Prior](#eq:p6)       | [Likelihood](#eq:p7) | [Posterior](#eq:p8) |
|---------------------------|-------------------------|------------------------|-----------------------|
| $\textnormal{Mode}(\theta)$ | $0.039$                 | $0.0047$               | $0.005$               |
| $\textnormal{E}[\theta]$              | $0.0476$                | --                     | $0.00517$             |
| $\textnormal{SD}(\theta)$ | $0.0207$                | --                     | $7.1\times 10^{-4}$   |
```


##### Example 2.7

 \
You are attempting to measure the **speed of a particle** $\theta$ and
quantify the uncertainty in the measurement. The experiment is set-up
such that the particle **travels** $1$km and the measurements are the
speed of the particle $X_i|\theta$, for $i = 1,\ldots, n$ in kilometers
per second ($km/s$).

There are a **multitude of factors** that additively result in
**measurement error** for the speed of the particle.

-   Using a **normal prior**, quantify **your prior beliefs** about
    $\theta$.

```{dropdown} Solution

Since there are physical constraints on the speed $\theta$, note that a
normal prior may not be appropriate. This is because a normal
distribution places probability mass over all the real line. Since we
are measuring speed, we should have $0 < \theta < c$, where
$c = 299792 km/s$ is the speed of light.

In any case, since we are asked to use a normal prior and we have no
knowledge of the type of the particle being considered, a prior with
little information is appropriate. For example:

$$
\theta \sim \mathcal{N}(1.5\times 10^5, (4\times 10^4)^2).$$

 This
places little mass outside the regions $\theta < 0$ and $\theta > c$ and
does not provide much information about the speed of the particle.

```

-   Why is a **normal distribution** appropriate for the likelihood?
    Write down the likelihood.

```{dropdown} Solution

Since there an accumulation of random errors, by the central limit
theorem $X_i|\theta\sim\mathcal{N}(\theta, \sigma^2)$ is a sensible
model. There will also be independence between each measurement. From
Example 1.7, the likelihood is

$$
f(\underline{x}|\theta,\sigma) =  (2\pi)^{-n/2}\sigma^{-n}
        \exp\left\{-\frac{1}{2\sigma^2}
        \left(\sum_{i=1}^n x_i^2-2\theta\sum_{i=1}^n x_i+n\theta^2\right)\right\}.$$

```

-   Assuming that the **standard deviation** of the measurements
    $\sigma$ is know, derive the posterior using a general
    $\theta \sim \mathcal{N}(\mu_0, \sigma^2_0)$ prior.

```{dropdown} Solution

We will use the fact that $e^{x+y} = e^x e^y$ throughout this example.
First, note that

$$
f(\underline{x}|\theta, \sigma) \propto \exp\left\{-\frac{1}{2\sigma^2}(-2\theta n\bar{x} + n\theta^2)\right\}$$

and

$$
\pi(\theta) \propto \exp\left\{-\frac{1}{2\sigma^2_0}(\theta^2 - 2 \mu_0 \theta)\right\}.$$

Therefore, 

$$
\begin{aligned}
        \pi(\theta|\underline{x}) &\propto \exp\left\{-\frac{1}{2\sigma^2}(-2\theta n\bar{x} + n\theta^2)\right\}\exp\left\{-\frac{1}{2\sigma^2_0}(\theta^2 - 2 \mu_0 \theta)\right\} \\
        &\propto \exp\left\{-\frac{1}{2}\left(\frac{n}{\sigma^2} \theta^2 - 2 \frac{n\bar{x}}{\sigma^2}\theta + \frac{\theta^2}{\sigma_0^2} - 2\frac{\mu_0}{\sigma_0^2} \theta \right)\right\} \\
        &\propto \exp\left\{-\frac{1}{2}\left(\left(\frac{n}{\sigma^2} + \frac{1}{\sigma_0^2}\right)\theta^2 - 2\left(\frac{n\bar{x}}{\sigma^2} + \frac{\mu_0}{\sigma_0^2}\right)\theta\right)\right\} \\
        &\propto \exp\left\{-\frac{1}{2}\left(\frac{\theta^2 - 2\left(\frac{n}{\sigma^2} + \frac{1}{\sigma_0^2}\right)^{-1}\left(\frac{n\bar{x}}{\sigma^2} + \frac{\mu_0}{\sigma_0^2}\right)\theta }{\left(\frac{n}{\sigma^2} + \frac{1}{\sigma_0^2}\right)^{-1}}\right)\right\}
    
\end{aligned}$$

 Let
$V = \left(\frac{n}{\sigma^2} + \frac{1}{\sigma_0^2}\right)^{-1}$ and
$M = \left(\frac{n}{\sigma^2} + \frac{1}{\sigma_0^2}\right)^{-1}\left(\frac{n\bar{x}}{\sigma^2} + \frac{\mu_0}{\sigma_0^2}\right)$.
Then, 

$$
\begin{aligned}
        \pi(\theta|\underline{x}) &\propto \exp\left\{-\frac{1}{2V}(\theta^2 - 2M\theta)\right\} \\
        &\propto \exp\left\{-\frac{1}{2V}(\theta^2 - 2M\theta)\right\} \exp\left(-\frac{1}{2V}M^2\right) \\
        &\propto \exp\left\{-\frac{1}{2V}(\theta^2 - 2M\theta + M^2)\right\} \\
        &\propto \exp\left\{-\frac{1}{2V}(\theta - M)^2\right\}.
    
\end{aligned}$$

 We recognise this as the form of a Normal distribution
with mean $M$ and variance $V$. So,

$$
\theta |\underline{x} \sim \mathcal{N}(M,V).$$

```

-   Suppose you attempted to measure the speed of the particle $4$ times
    and obtained the measurements (in $km/s$): 306135, 293227, 307985,
    301298.

    Using your prior, setting $\sigma = 5000$ and using this data,
    derive your **posterior distribution**.

```{dropdown} Solution

Using our posterior update rule from part (c) with
$\bar{x} = 302161.25$, $\sigma^2 = 5000^2$, $\mu_0 = 150000$ and
$\sigma_0^2 = 40000^2$, the posterior is

$$
\theta | \underline{x} \sim \mathcal{N}(271729, 2236.068^2).$$

```

**Summary:**

From part (c), we know that if our data $\underline{x}$ is IID normally
distributed with unknown mean $\theta$ and we choose a Normal prior
distribution 

$$
\theta \sim \mathcal{N}(\mu_0,\sigma^2_0)$$

 to quantify
our uncertainty about $\theta$, then the distribution of the posterior
$\theta|\underline{x}\sim \mathcal{N}(M, V)$ is again normally
distributed with mean $M$ and variance $V$, where 

$$
\begin{aligned}
    M &= \left(\frac{n}{\sigma^2} + \frac{1}{\sigma_0^2}\right)^{-1}\left(\frac{n\bar{x}}{\sigma^2} + \frac{\mu_0}{\sigma_0^2}\right), \\
    V &= \left(\frac{n}{\sigma^2} + \frac{1}{\sigma_0^2}\right)^{-1}.
\end{aligned}$$

The prior, likelihood function and posterior can be seen in Figure .

![Prior, likelihood function and posterior of the speed of the particle
$\theta$.](images/speed_priorposterior.svg)

##### Example 2.8

 \
Let $Y$ be the retreat, in feet, of the *Zachariae Isstrøm* glacier. A
*Pareto* distribution with rate $\theta$ is often used to model such
geophysical activity, with probability density function

$$
f(y|\kappa,\theta) = \theta\kappa^{\theta}y^{-(\theta+1)}, \qquad \theta,\kappa>0 \mathrm{ and } y>\kappa.$$

-   Obtain the likelihood function for $\theta$ given the parameter
    $\kappa$ and some observed data $y_{1}, y_{2}, \ldots, y_{n}$
    (independent).

```{dropdown} Solution

The likelihood function is simply the product of the probability density
function evaluated at each observation $y_{i}$, ($i=1, \ldots, n$), i.e.

$$
\begin{aligned}
    f(\underline{y} | \theta, \kappa) &= \theta\kappa^{\theta}y_{1}^{-(\theta+1)} \times \cdots \times \theta\kappa^{\theta}y_{n}^{-(\theta+1)}  \notag \\
                                   & \notag \\
                                   &= \theta^{n} \kappa^{n\theta} \prod_{i=1}^{n}y_{i}^{-(\theta+1)}.
    
\end{aligned}$$

```

-   Suppose we observe a retreat of 20 feet at the *Zachariae Isstrøm*
    glacier in 2012. Write down the likelihood function for $\theta$.

```{dropdown} Solution

We simply substitute $n=1$ and $y_{1}=20$ into the likelihood, giving

$$
f(\theta|\kappa,y_{1}=20)     =\theta \kappa^{\theta} 20^{-(\theta+1)}.$$

```

-   Using the prior $\theta \sim \mathrm{Gamma}(9,0.36)$ for the rate of
    retreat and assuming $\kappa$ is known to be 12, obtain the
    posterior distribution $\pi(\theta|y_{1}=20)$.

```{dropdown} Solution

Using Bayes Theorem, and following the examples in Chapter 2, we know
that

$$
\pi(\theta|y_{1}=20) \propto \pi(\theta)\times f(y_1 =20 | \theta).$$

Recall from Example 3.4 that our elicited prior for $\theta$ is
$\mathrm{Gamma}(9,0.36)$, which has density

$$
\pi(\theta) = \frac{0.36^{9}\theta^{8}e^{-0.36\theta}}{\Gamma(9)}.$$

Combining this with the likelihood above (and using $\kappa=12$) gives

$$
\begin{aligned}
        \pi(\theta|y_{1}=20)   &= \frac{0.36^{9}\theta^{8}e^{-0.36\theta}}{\Gamma(9)} \times \theta 12^{\theta}20^{-(\theta+1)} \notag \\
           &  \notag \\
                               &\propto \theta^{9} e^{-0.36 \theta} 12^{\theta} 20^{-(\theta+1)}  \notag \\
            & &\notag   \\
                               &\propto \theta^{9}e^{-0.36 \theta}12^{\theta} 20^{-\theta}.
        
\end{aligned}$$

Now consider the term $12^{\theta}20^{-\theta}$. Taking logs, we get

$$
\theta \mathrm{ln}12 - \theta \mathrm{ln}20  = (\mathrm{ln}12 - \mathrm{ln}20)\theta;$$

exponentiating to 're--balance', you should see that

$$
12^{\theta}20^{-\theta} = e^{(\mathrm{ln}12 - \mathrm{ln}20)\theta}.$$

Substituting back into (3.6) gives 

$$
\begin{aligned}
          \pi(\theta|y_{1}=20)  &\propto& \theta^{9}e^{-0.36 \theta} e^{(\mathrm{ln}12 - \mathrm{ln}20)\theta}     \qquad \mathrm{i.e.}     \\
                   & & \\
                                &\propto& \theta^{9}e^{-0.36\theta +(\mathrm{ln}12-\mathrm{ln}20)\theta}    \\
                         & & \\
                                &\propto& \theta^{9}e^{-0.87 \theta}.
        
\end{aligned}$$

 Referring to our definition of the gamma distribution on
page 25 of these notes, you should recognise this as a
$\mathrm{Gamma}(10,0.87)$ distribution.

```

![Prior (dashed) and posterior (solid) densities for the rate of glacial
retreat at the *Zachariae Isstrøm*
glacier.](images/glacier_posterior.svg)

## Section 2.3: Conjugacy

In Example 2.5, a gamma prior leads to a gamma posterior. This is an
example of *conjugacy*, where choosing a prior in a family of
distributions always leads to a posterior in the same family. The formal
definition is as follows:

##### Definition 2.4: Conjugate Prior

 \
Suppose that data $\underline{x}$ are modelled with distribution
$f(\underline{x}|\theta)$. A family $\mathcal{P}$ of prior distributions
for $\theta$ is *conjugate* to $f(\underline{x}|\theta)$ if for every
prior distribution $\pi(\theta)\in\mathcal{P}$, the posterior
distribution $\pi(\theta|\underline{x})$ is also in $\mathcal{P}$.

Notice that the conjugate family depends crucially on the model chosen
for the data $\underline{x}$. So we say (for example) that the Gamma
distribution is a conjugate prior for the exponential model.

There are usually simple formulae to update the conjugate prior to the
corresponding posterior. For instance, Example 2.5 showed that a
$\mathrm{Gamma}(g,h)$ prior for an exponential model leads to a
$\mathrm{Gamma}(g+n,h+n\bar x)$ distribution. For this reason, conjugate
priors are usually a convenient choice to reduce the mathematical and
computational effort of deriving the posterior.

Some examples of conjugacy are as follows:

1.  Binomial random sample, Beta prior distribution $\longrightarrow$
    Beta posterior distribution (Examples 2.1 -- 2.3)

2.  Exponential random sample, Gamma prior distribution
    $\longrightarrow$ Gamma posterior distribution (Examples 2.4 -- 2.5)

3.  Normal random sample (known variance), Normal prior distribution
    $\longrightarrow$ Normal posterior distribution (Example 2.6)

4.  $\mathrm{Gamma}(k,\theta)$ random sample ($k$ known), Gamma prior
    distribution $\longrightarrow$ Gamma posterior distribution.

##### Example 2.9

 \
Suppose we have a random sample from a normal distribution. In Bayesian
statistics, when dealing with the normal distribution, the mathematics
is more straightforward if we work with the precision
($=1/\mathrm{variance}$) of the distribution rather than the variance
itself. So we will assume that this population has unknown mean $\mu$
but known precision $\tau$. That is,

$$
X_i|\mu\sim \mathcal{N}(\mu,1/\tau),$$

 for $i=1,2,\ldots,n$
(independent), where $\tau$ is known. Suppose our prior beliefs about
$\mu$ can be summarised by a $\mathcal{N}(b,1/d)$ distribution, with
probability density function 

```{math}
:label: eq:p11
\pi(\mu)=
\left(\frac{d}{2\pi}\right)^{1/2}\exp\left\{-\frac{d}{2}(\mu-b)^2\right\}.
```

Determine the posterior distribution for $\mu$.

```{dropdown} Solution

The likelihood function for $\mu$ is 

$$
\begin{aligned}
        f(\underline{x}|\mu)&=\prod_{i=1}^n 
        \left(\frac{\tau}{2\pi}\right)^{1/2}
        \exp\left\{-\frac{\tau}{2}(x_i-\mu)^2\right\}  \\
        &= \left(\frac{\tau}{2\pi}\right)^{n/2}
        \exp\left\{-\frac{\tau}{2}\sum_{i=1}^n (x_i-\mu)^2\right\} \\
        &= \left(\frac{\tau}{2\pi}\right)^{n/2}
        \exp\left\{-\frac{\tau}{2}\sum_{i=1}^n (x_i-\bar x+\bar x-\mu)^2\right\} \\
        &= \left(\frac{\tau}{2\pi}\right)^{n/2}
        \exp\left\{-\frac{\tau}{2}\left[\sum_{i=1}^n (x_i-\bar x)^2+n(\bar x-\mu)^2\right]\right\} 
        
\end{aligned}$$

 Let 

$$
s^2 =\frac{1}{n}\sum_{i=1}^n (x_i-\bar x)^2$$

 and
then 

<div style="text-align: center;">
    <a id="eq:p10"></a>
    <div class="math notranslate nohighlight">
        $$
f(\underline{x} | \mu) = \left(\frac{\tau}{2\pi}\right)^{n/2}
        \exp\left\{-\frac{n\tau}{2}\left[s^2+(\bar x-\mu)^2\right]\right\}.
        $$
    </div>
</div>

Applying Bayes Theorem, the posterior density function is

$$
\begin{aligned}
        \pi(\mu|\underline{x})&\propto\pi(\mu)f(\underline{x} | \mu)\\
        &\propto \left(\frac{d}{2\pi}\right)^{1/2}\exp\left\{-\frac{d}{2}(\mu-b)^2\right\} \times \left(\frac{\tau}{2\pi}\right)^{n/2}
        \exp\left\{-\frac{n\tau}{2}\left[s^2+(\bar x-\mu)^2\right]\right\} \\
        &=k_1\exp\left\{-\frac{1}{2}\left[d(\mu-b)^2+n\tau(\bar x-\mu)^2\right]\right\}
        
\end{aligned}$$

 where $k_1$ is a constant that does not depend on $\mu$.
Now the exponent can be simplified by expanding terms in $\mu$ and then
completing the square, as follows.

We have 

$$
\begin{aligned}
d(\mu-b)^2&+n\tau(\bar x-\mu)^2\\
&=d(\mu^2-2b\mu+b^2)+n\tau(\bar x^2-2\bar x\mu+\mu^2)\\
&=(d+n\tau)\mu^2-2(db+n\tau\bar x)\mu
+db^2+n\tau\bar x^2\\
&=(d+n\tau)\left\{\mu
-\left(\frac{db+n\tau\bar x}{d+n\tau}\right)\right\}^2+c
\end{aligned}$$

 where $c$ does not depend on $\mu$. Let

<div style="text-align: center;">
    <a id="eq:p12"></a>
    <div class="math notranslate nohighlight">
        $$
B=\frac{db+n\tau\bar x}{d+n\tau}\quad\quad\text{and}\quad\quad
D=d+n\tau.
        $$
    </div>
</div>

 Then 

<div style="text-align: center;">
    <a id="eq:p13"></a>
    <div class="math notranslate nohighlight">
        $$
\begin{aligned}
\pi(\mu|\underline{x})&=k_1\exp\left\{-\frac{D}{2}(\mu-B)^2-\frac{c}{2}\right\}
\notag \\
&=k\exp\left\{-\frac{D}{2}(\mu-B)^2\right\},
\end{aligned}
        $$
    </div>
</div>

 where $k$ is a constant that does not depend on $\mu$.
Therefore, the posterior distribution takes the form
$k\exp\{-D(\mu-B)^2/2\}$, $-\infty<\mu<\infty$ and so must be a normal
distribution: we have $\mu|\underline{x}\sim \mathcal{N}(B,1/D)$.

```

**Summary:**

If we have a random sample from a $\mathcal{N}(\mu,1/\tau)$ distribution
(with $\tau$ known) and our prior beliefs about $\mu$ follow a
$\mathcal{N}(b,1/d)$ distribution then, after incorporating the data,
our (posterior) beliefs about $\mu$ follow a $\mathcal{N}(B,1/D)$
distribution.

Notice that the way prior information and observed data combine is
through the parameters of the normal distribution:

$$
B = \frac{db+n\tau\bar x}{d+n\tau}
\quad\quad\text{and}\quad\quad 
D = d+n\tau.$$

 Notice also that the posterior variance (and precision)
does not depend on the data, and the posterior mean is a convex
combination of the prior and sample means, that is,

$$
B=\alpha b+(1-\alpha)\bar x,$$

 for some $\alpha\in(0,1)$. This
equation for the posterior mean, which can be rewritten as

$$
\text{E}[\mu|\underline{x}]=\alpha \text{E}[\mu]+(1-\alpha)\bar x,$$

arises in other models and is known as the *Bayes linear rule*.

The changes in our beliefs about $\mu$ are summarised in Table . Notice
that the posterior mean is greater than the prior mean if and only if
the likelihood mode (sample mean) is greater than the prior mean, that
is

$$
\text{E}[\mu|\underline{x}] > \text{E}[\mu] \quad\iff\quad \textnormal{Mode}(f(\underline{x}|\mu))>\text{E}[\mu].$$

Also, the standard deviation of the posterior distribution is smaller
than that of the prior distribution.

```{table} Changes in beliefs about $\mu$.
:widths: auto
:align: center
:name: tab:norknown

|                            | [Prior](#eq:p11)       | [Likelihood](#eq:p10) | [Posterior](#eq:p13) |
|----------------------------|--------------------------|-------------------------|------------------------|
| $\textnormal{Mode}(\mu)$   | $b$                      | $\bar x$                | $(db+n\tau\bar x)/(d+n\tau)$ |
| $\textnormal{E}[\mu]$                  | $b$                      | --                      | $(db+n\tau\bar x)/(d+n\tau)$ |
| $\textnormal{Precision}(\mu)$ | $d$                   | --                      | $d+n\tau$              |
```


##### Example 2.10

 \
The ages of *Ennerdale granophyre* rocks can be determined using the
relative proportions of rubidium--87 and strontium--87 in the rock. An
expert in the field suggests that the ages of such rocks (in millions of
years) $X|\mu\sim \mathcal{N}(\mu,8^2)$ and that a prior distribution
$\mu\sim \mathcal{N}(370,20^2)$ is appropriate. A rock is found whose
chemical analysis yields $x=421$. What is the posterior distribution for
$\mu$ and what is the probability that the rock will be older than 400
million years?

```{dropdown} Solution

We have $n=1$, $\bar x=x=421$, $\tau=1/64$, $b=370$ and $d=1/400$.
Therefore, we have: 

$$
\begin{aligned}
        B&=\frac{db+n\tau\bar x}{d+n\tau}
        =\frac{370/400+421/64}{1/400+1/64}=414.0\\
        D&=d+n\tau=1/400+1/64=1/7.43^2
        
\end{aligned}$$

and so the posterior distribution is
$\mu|x=421\sim \mathcal{N}(414.0,7.43^2)$. The (posterior) probability
that the rock will be older than 400 million years is 

$$
\begin{aligned}
        Pr(\mu>400|x=421)=0.9702
        
\end{aligned}$$

calculated using the $R$ commands $1-\texttt{pnorm}(400,414,7.43)$ or
$1-\texttt{pnorm}(-1.884)$ or $\texttt{pnorm}(1.884)$. Without the
chemical analysis, the only basis for determining the age of the rock is
via the prior distribution: the (prior) probability that the rock will
be older than 400 million years is $Pr(\mu>400)=0.0668$ calculated using
the R command $1-\texttt{pnorm}(400,370,20)$.

```

This highlights the benefit of taking the chemical measurements. Note
that the large difference between these probabilities is not necessarily
due to the expert's prior distribution being inaccurate, *per se*, it is
probably due to the large prior uncertainty about rock ages, as shown in
Figure .

![Prior (dashed) and posterior (solid) densities for the age of the
rock.](images/priorposterior4.svg)
