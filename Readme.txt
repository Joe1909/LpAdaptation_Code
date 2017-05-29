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

the data used in the paper can be found in the folder DataPaper

Figure3: 
	out_Storn10.mat
	Plot_Figure3.m
		 
Figure4: 
	out_Folium.mat, out_Handle.mat

Figure5: 
	all out_Lp_stretched20D_*.mat files
	Plot_Figure5.m
		
Figure6: 
	out_SCfilter.mat
	Plot_Figure6.m
		 
two component system: 
	all out_TCS_*.mat files	 		 		 

