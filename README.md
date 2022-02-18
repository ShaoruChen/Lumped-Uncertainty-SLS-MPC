# Lumped-Uncertainty-SLS-MPC

This repo. contains the codes for implementing the robust model predictive contro (MPC) methods used in the paper [System Level Synthesis-based Robust Model Predictive
Control through Convex Inner Approximation](https://arxiv.org/pdf/2111.05509.pdf). It aims to solve robust MPC of linear time-invariant systems subject to both norm-bounded multiplicative model uncertainty and additive disturbances:

![](https://latex.codecogs.com/svg.image?x_{t&plus;1}&space;=&space;(\hat{A}&space;&plus;&space;\Delta_A)x_t&space;&plus;&space;(\hat{B}&space;&plus;&space;\Delta_B)&space;u_t&space;&plus;&space;w_t)  

with ![](https://latex.codecogs.com/svg.image?\inline&space;\lVert&space;\Delta_A&space;\rVert_\infty&space;\leq&space;\epsilon_A,&space;\lVert&space;\Delta_B&space;\rVert_\infty&space;\leq&space;\epsilon_B,&space;\lVert&space;w_t&space;\rVert_\infty&space;\leq&space;\sigma_w).

## Summary
At each time instnat in MPC, the proposed method uses [system level synthesis](https://arxiv.org/abs/1904.01634) (SLS) to solve the robust optimal control problem in the space of closed-loop **system responses** which allows a novel constraint tightening procedure that would otherwise be impossible to apply. The proposed method, which we denote as lumped uncertainty SLS MPC, demonstrates significant improvement in conservatism compared with other competitive robust MPC baselines across a wide range of numerical examples.

## Content
- This repo. contains implementation of the proposed method, [lumped uncertainty SLS MPC](https://github.com/ShaoruChen/Lumped-Uncertainty-SLS-MPC/blob/1fd27cbc27ee93a9ef0f1cb68f9c21000298ed2e/mpc/SLSMPC.m#L1015), together with three robust MPC baselines, namely [tube-MPC](https://github.com/ShaoruChen/Lumped-Uncertainty-SLS-MPC/blob/1fd27cbc27ee93a9ef0f1cb68f9c21000298ed2e/mpc/SLSMPC.m#L655) [1], [disturbance feedback MPC using uniform model uncertainty abstraction](https://github.com/ShaoruChen/Lumped-Uncertainty-SLS-MPC/blob/1fd27cbc27ee93a9ef0f1cb68f9c21000298ed2e/mpc/SLSMPC.m#L869) [2], and [SLS MPC using grid search of hyperparameters](https://github.com/ShaoruChen/Lumped-Uncertainty-SLS-MPC/blob/1fd27cbc27ee93a9ef0f1cb68f9c21000298ed2e/mpc/SLSMPC.m#L611) [3] . All these methods can handle both model uncertainty and additive disturbances. 

    [1] Langson, Wilbur, Ioannis Chryssochoos, S. V. Raković, and David Q. Mayne. [Robust model predictive control using tubes.](https://www.sciencedirect.com/science/article/abs/pii/S0005109803002838) Automatica 40, no. 1 (2004): 125-133.

    [2] Bujarbaruah, Monimoy, Ugo Rosolia, Yvonne R. Stürz, and Francesco Borrelli. [A simple robust MPC for linear systems with parametric and additive uncertainty.](https://arxiv.org/abs/2103.12351) In 2021 American Control Conference (ACC), pp. 2108-2113. IEEE, 2021.

    [3] Chen, Shaoru, Han Wang, Manfred Morari, Victor M. Preciado, and Nikolai Matni. [Robust closed-loop model predictive control via system level synthesis.](https://arxiv.org/abs/1911.06842) In 2020 59th IEEE Conference on Decision and Control (CDC), pp. 2152-2159. IEEE, 2020.

- Numerical examples in the paper [System Level Synthesis-based Robust Model Predictive
Control through Convex Inner Approximation](https://arxiv.org/pdf/2111.05509.pdf) can be produced by running the codes under the [example](https://github.com/ShaoruChen/Lumped-Uncertainty-SLS-MPC/tree/main/examples) folder. 

## Installation
Add the [mpc](https://github.com/ShaoruChen/Lumped-Uncertainty-SLS-MPC/tree/main/mpc) folder to MATLAB path and then you can run the [examples](https://github.com/ShaoruChen/Lumped-Uncertainty-SLS-MPC/tree/main/examples) in the paper. 

### Required toolboxes
[Yalmip](https://yalmip.github.io/) for formulating the control problems. [MOSEK](https://docs.mosek.com/9.3/toolbox/install-interface.html) is used as the default solver in the codes. 

[MPT3](https://www.mpt3.org/) for polyhedron operations. 

[MatlabProgressBar](https://www.mathworks.com/matlabcentral/fileexchange/57895-matlabprogressbar) for progress display (Not required if you remove the progress function in each for-loop, e.g. for i = progress(1:10) --> for i = 1:10).
