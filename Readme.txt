Code for design centering and volume estimation as used and described in

J. Asmus, C. L. Mueller and I. F. Sbalzarini. Lp-Adaptation: Simultaneous 
Design Centering and Robustness Estimation of Electronic and Biological
Systems. Scientific Reports 2017

LpAdaptation needs 2 inputs: 
1) an oracle, describing the specifications that need to be met. 
   The oracle is a binary function, yielding 1, if the specifications are met and 0 if they are not
2) a starting point for which the oracle yields 1

in folder oracles, there are different examples of oracles

additionally to these two inputs, different options can be passed, see LpAdaptation/getDefaultOptions.m 
those options are given in a struct, LpAdaptation is then called as
	LpAdaptation(oracle,xstart,inopts)
see also examples in folder examples

the data used in the paper can be found here: 
	https://asmus@git.mpi-cbg.de/asmus/LpAdaptation_Code.git
