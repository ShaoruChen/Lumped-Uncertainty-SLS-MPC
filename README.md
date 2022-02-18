# Lumped-Uncertainty-SLS-MPC

This repo. contains the codes for implementing the robust model predictive contro (MPC) methods used in the paper [System Level Synthesis-based Robust Model Predictive
Control through Convex Inner Approximation](https://arxiv.org/pdf/2111.05509.pdf). It aims to solve robust MPC of linear time-invariant systems subject to both norm-bounded multiplicative model uncertainty and additive disturbances:

![](https://latex.codecogs.com/svg.image?x_{t&plus;1}&space;=&space;(\hat{A}&space;&plus;&space;\Delta_A)x_t&space;&plus;&space;(\hat{B}&space;&plus;&space;\Delta_B)&space;u_t&space;&plus;&space;w_t)  

with ![](https://latex.codecogs.com/svg.image?\inline&space;\lVert&space;\Delta_A&space;\rVert_\infty&space;\leq&space;\epsilon_A,&space;\lVert&space;\Delta_B&space;\rVert_\infty&space;\leq&space;\epsilon_B,&space;\lVert&space;w_t&space;\rVert_\infty&space;\leq&space;\sigma_w).

## Required toolboxes
[Yalmip](https://yalmip.github.io/) 

[MPT3](https://www.mpt3.org/)

[MatlabProgressBar](https://www.mathworks.com/matlabcentral/fileexchange/57895-matlabprogressbar) (Not required if you remove the progress function in each for-loop, e.g. for i = progress(1:10) --> for i = 1:10).
