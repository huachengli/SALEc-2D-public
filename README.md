#### SALEc2d

SALEc2d is a 2d version of [SALEc](https://github.com/huachengli/SALEc-public) aimed to calculate the process of impact craters with higher resolution.

Features implemented in SALEc2d:

1. Parallelism/SALEc2d can used MPI to accelerated the simulation.

2. Absorbing Boundary/Sponge boundary condition(SBC) and Perfect Match layers(PML) boundary have been implemented. The layers recommended for SBC and PML is 50 and 30 respectively.

3. A mixed strategy for constructing material internal surface based on volume of fluid(VOF). The distance between vertexes and vacuum-materials is calculated.

4. Improved Stability/Artificial diffusion is applied in restricted section. Materials isolated which occupy less than one element has been deleted. Some void vertexes get velocity from its neighbor elements based on its distance from the vacuum-materials interface.

see: 

 	1. [Parallel numerical simulation of impact crater with perfect matched layers](https://arxiv.org/abs/2403.04267)
 	2. [Analysis on the source position of Zhinyu crater ejecta](https://doi.org/10.1016/j.icarus.2025.116579)

Author: 
    Li Huacheng, huacheng_li@pku.edu.cn
