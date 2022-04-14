# coeff-var
Code for computing coefficient of variation for small samples

This code computes the coefficient of variation (CV) and some other stats for small samples (indicated by the * added to CV) 
for a given set of measurements which are assumed to be for the same or similar object, using the same measurand. 
Stats are adjusted for small sample size. Paper ref: Belz, Popovic & Mille (2022) Quantified Reproducibility Assessment of NLP Results,
ACL'22.

In this self-contained version, the set of measurements on which CV is computed is assigned to the variable set_of_set_of_measurements
(see examples in code below).

The reproducibility stats reported in the output are: 
* the unbiased coefficient of variation
* the sample mean
* the unbiased sample standard deviation with 95% confidence intervals, estimated on the basis of the standard error of the unbiassed sample variance
* the sample size
* the percentage of measured valued within two standard deviations
* the percentage of measured valued within one standard deviation

Example narrative output:

The unbiased coefficient of variation is 1.5616560359100269 \
for a mean of 85.58285714285714 , \
unbiased sample standard deviation of 1.2904233075765223 with 95\% CI (0.4514829817654973, 2.1293636333875474) ,\
and a sample size of 7 . \
100.0 % of measured values fall within two standard deviations. \
71.429 % of measured values fall within one standard deviation. 

NOTE:
* CV assumes all measurements are positive; if they're not, shift measurement scale to start at 0
* for fair comparison across studies, measurements on a scale that doesn't start at 0 need to be shifted to a scale that does start at 0 
