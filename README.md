#### SALEc2d

SALEc2d is a 2d version of [SALEc](https://github.com/huachengli/SALEc-public) aimed to calculate the process of impact craters with higher resolution.

Features implemented in SALEc2d:

1. Parallelism/SALEc2d can used MPI to accelerated the simulation.

2. Absorbing Boundary/Sponge boundary condition(SBC) and Perfect Match layers(PML) boundary have been implemented. The layers recommended for SBC and PML are 50 and 30 respectively.

3. Improved Stability.

see: 

[Parallel numerical simulation of impact crater with perfect matched layers](https://arxiv.org/abs/2403.04267)

[Analysis on the source position of Zhinyu crater ejecta](https://doi.org/10.1016/j.icarus.2025.116579)

Author: 
    Li Huacheng, huacheng_li@pku.edu.cn
