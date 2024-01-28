# MAS2903: Introduction to Bayesian Methods

**Dr. Matthew Fisher**

Office: Room 4.13, Herschel Building

Email: `matthew.fisher@newcastle.ac.uk`

**Assessment**

MAS2903 is a 10 credit module assessed via a formal written examination
in the summer (likely a traditional 2 hour written paper, taken
in-person under exam conditions) and some in-course assessment. The
make-up is 80% summer exam, 20% in-course assessment. Further details:

-   Summer exam, 100 marks in total

    -   Section A (40 marks): Consisting mainly of short,
        straightforward, data-response style questions; perhaps with
        some multiple choice elements

    -   Section B (60 marks): Longer, harder questions, requiring some
        contextual discussion and interpretation. Questions more
        mathematically challenging

    -   Answer *all* questions in *both* parts

-   In-course assessment, 80 marks in total

    -   Problem Solving Exercises 1 (20 marks)

    -   Numbas Assessment 1 (20 marks)

    -   Problem Solving Exercises 2 (20 marks)

    -   Numbas Assessment 2 (20 marks)

Some of the questions in the in-course assessment will require you to
produce plots in `R` or `Python`, so it might be a good idea to type up your
solutions and upload a single PDF file (which includes any necessary
plots/computing work) to Canvas.

**Books and other resources**

At various points in the module you will be directed to additional
reading, should you want extra examples or just another source of
information relating to the material presented in the notes. This will
primarily be text books (many of which are available to view online),
but perhaps other sources (e.g. journal publications, magazine articles
etc.).

Some useful books include:

-   "*Bayes Rule: A Tutorial Introduction to Bayesian Analysis*" - James
    Stone

-   "*Doing Bayesian Data Analysis: A Tutorial with R, JAGS, and
    Stan*" - John Krushke

-   "*Bayesian Statistics: An Introduction*" - Peter Lee

"Bayes rule" is a good introduction to the main concepts in Bayesian
statistics but doesn't cover everything in this course. The other books
are broader references which go well beyond the contents of this course.

**Acknowledgements**

These notes are a small modification on notes by Dr Lee Fawcett.

# Preface: What is Bayesian Statistics?

There are two main approaches to statistics: *frequentist* (or
classical) statistics and *Bayesian* statistics. All of the statistics
teaching you've encountered so far is likely to be about frequentist
methods. Bayesian methods are substantially different and can feel quite
strange to start with. So, before starting the main course material,
this *non-examinable* preface explains the concepts behind Bayesian
statistics and how they differ from frequentist approaches.

One way of defining statistics is as a way to learn about the world from
some data which is subject to random variation. For example, in *climate
science* we want to learn about the climate given some imperfect
measurements. Some statistical questions that we might ask in this field
are:

-   What is our best estimate of the world temperature this year?\
    This is referred to as a *point estimation* problem.

-   What is a plausible range of values?\
    This is referred to as an *interval estimation* (or uncertainty
    quantification) problem.

-   What will the temperature be in 100 years?\
    This is a *prediction* problem.

-   Is the climate warming?\
    This is a *hypothesis testing* problem.

Frequentist and Bayesian methods address all these types of problem:
point estimation, interval estimation, prediction and hypothesis
testing. But they use different approaches to do so. Some typical
frequentist approaches include:

-   Least squares (point estimation)

-   Maximum likelihood (point estimation)

-   Confidence intervals (interval estimation)

-   Test statistics and p-values (hypothesis testing)

Bayesian statistics doesn't use any of these familiar methods! However,
there are some connections between them and their Bayesian alternatives
which will be explored later in the course. For example we'll see that
the idea of the *likelihood function* is very important in Bayesian
statistics.

## Motivating example

Air France Flight 447 disappeared over the Atlantic Ocean on June 1st
2009, during an overnight flight from Rio de Janeiro to Paris. A search
and rescue mission was launched the following morning. After five days
some floating wreckage was found, but there was no sign of the "black
box" recorders which would reveal what had happened to the flight, and
the search stalled.

The black boxes could be anywhere within an area of the South Atlantic
the size of Switzerland (6,500 square miles). The "mid-Atlantic ridge"
on the ocean floor lies between two tectonic plates and is just as
mountainous as Switzerland! It's also so remote that scientists have not
yet charted the sea-bed. The search method used was sonar-detectors,
emitting sound waves which would bounce back once they hit something.
However, after two years of meticulously searching an area north of the
plane's trajectory (after analysing debris drift), nothing had been
found.

Metron, Inc., of Retson, Virginia had been performing Bayesian analyses
to help the search effort. Included in their analysis were:

1.  Data from 9 previous airline crashes involving loss of pilot
    control - reduced search area to 1,600 square miles.

2.  Expert opinions on the credibility of the flight data.

3.  Expert opinions about whether or not the black box 'pingers' might
    have become damaged on impact.

4.  Positions/recovery times of bodies found drifting - expert opinions
    assigned to the reliability of this data because of the turbulent
    equatorial waters.

5.  Expert information from oceanographers: Sea state, visibility,
    underwater geography,...

All were combined through Bayes Theorem to give a probability map of the
most likely locations to search. In April 2011 they decided to try
assuming that the 'pingers' in item 3 were probably damaged. One week
later the black boxes were found![^1]

## The frequentist vs. Bayesian debate

This example illustrates some key features of the Bayesian approach:

-   It incorporates expert opinion.

-   It outputs a probability distribution over the possible choices.

-   The calculations required are all based on Bayes Theorem.

In contrast, in frequentist statistics, it is difficult (but not
impossible) to incorporate opinions, and the output is not in the form
of a probability distribution. For example a 95% frequentist confidence
interval is a range which will contain the true value 95% of the time if
the analysis is repeated a large number of times for different data
obtained under the same process[^2]. In contrast a Bayesian 95% interval
estimate (see Chapter 4) is a range with a 95% probability of containing
the true value. This Bayesian version is arguably much easier to
interpret.

There has been a decades-long argument about whether Bayesian or
frequentist methods are better. Both sides have some good points. For
example, there is the argument we've just made about Bayesian interval
estimates being easy to interpret. On the other hand one of the main
arguments against the Bayesian approach is that it is *subjective*. This
is because its results are based on expert opinions, and if we got a
different opinion from another expert our results would change even
though the data are the same! Bayesians have counter-arguments.

In practice many problems are better suited to be tackled with methods
from one approach or the other. For example Bayesian statistics was well
suited to the Air France example as incorporating expert knowledge was
particularly important in the solution. Therefore it's good to be
familiar with both Bayesian and frequentist methods.

## The history of Bayesian statistics

The reason why we use the word "Bayesian" is because Bayes Theorem is
crucial for statistical analysis if we adopt the Bayesian approach.
Thomas Bayes (1702--1761; see picture below) was a Presbyterian minister
in Tunbridge Wells, Kent. Bayes solution to a problem in probability
theory was presented in the *Essay towards Solving a Problem in the
Doctrine of Chances*, published posthumously by his friend Richard Price
in the *Philosophical Transactions of the Royal Society of London*. This
paper gave us Bayes Theorem. Its modern form is due to Pierre-Simon
Laplace who generalised Bayes work and produced the first theory of
Bayesian inference, which was referred to as "inverse probability".
Laplace's theory had weaknesses[^3] and in the early 20th century it was
replaced by the more rigorous frequentist approach. R. A. Fisher was one
of the leading creators of this approach. He was also the first to use
the term "Bayesian statistics".

Nonetheless there was still some interest in Bayesian methods. For
example, Alan Turing's group invented some new Bayesian methods to help
break the Enigma code in World War II. Before and, increasingly, after
the war a few researchers worked on developing a rigorous theory of
Bayesian statistics.

Bruno de Finetti (1906--1985) was an Italian probabilist who developed
his ideas on *subjective probability*[^4] from the 1920s onwards,
drawing upon ideas from H. Jeffreys, I.J. Good (one of the Enigma
codebreakers) and B.O. Koopman. His classic book on the topic is *The
Theory of Probability* (1974).

Dennis Lindley (1923--2013) was a leading British Bayesian. In his early
career, he worked to find a mathematical basis for the subject of
statistics. In 1954, Lindley met Leonard Savage and both found a deeper
justification for statistics in Bayesian theory, turning into critics of
the classical statistical inference they had hoped to justify.

Many of these theoretical advances are covered in this course. However,
Bayesian statistics was still little used in practice. The difficulty
was that applying Bayes theorem often became very mathematically
challenging, resulting in problems that were too difficult or impossible
to solve by hand.

### Recent developments

In the 1990s a solution to these mathematical problems was developed,
making use of the emergence of more powerful computers and "Markov chain
Monte Carlo" (MCMC) numerical algorithms. If you take 3rd and 4th year
Bayesian courses you will be introduced to these techniques. MCMC, and
subsequent developments, have revolutionised the use of Bayesian
Statistics, to the extent that today Bayesian data analyses are as
popular as frequentist approaches and are routinely used in fields as
diverse as artificial intelligence, biology, astrophysics and sociology.

### Newcastle's contribution

Since 1980, the number of academic staff in Mathematics & Statistics at
Newcastle publishing advanced research using Bayesian methods has
increased dramatically. In the 1980s, there was only one Bayesian at
Newcastle. Now there are at least 12. Our Bayesian research topics
include environmental extremes, genetics, and modelling biological
systems at a molecular level.

[^1]: For more details see
    `https://www.technologyreview.com/s/527506/how-statisticians-found-`
    `air-france-flight-447-two-years-after-it-crashed-into-atlantic/`.

[^2]: Here 95% is the *long-run frequency* that the interval contains
    the correct value. This kind of property is where the name
    "frequentist" comes from.

[^3]: For example it relied mainly on *flat priors*, which will be
    discussed and criticised in Chapter 3.

[^4]: Covered in Chapter 1.
