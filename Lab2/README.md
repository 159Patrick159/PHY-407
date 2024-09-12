# Lab 2: Numerical Errors and Integration
This sub-directory contains the python scripts and written report for Lab 2 for Computational Physics. The breakdown of the lab is the following

## Q1 - Understanding Numerical Errors from Standard Deviation Calculations
We examine the contraints of perfoming numerical operations on big numbers and how this can lead to inaccuracies in computations. We examine the effects by computing the standard deviation of two gaussians one with high mean and one with low mean. Then we use a two-pass standard deviation equation and a one-pass standard deviation equation. We find the two-pass equation to perform more accurately than the one-pass algorithm. Finally we suggest an alternative implementation for the one-pass algorithm to improve its accuracy.

<p align='center'>
    <img src="Figures/Q1cPlot.png" title="Generated histrograms with different means but same standard deviations." height="80%" width="80%">
</p>
