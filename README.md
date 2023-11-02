# EcoLand
EcoLand is a pipeline to simulate the stochastic dynamics of the ecological systems (Jinchao Lv, Jin Wang* and Chunhe Li*, Landscape quantifies the intermediate state and transition dynamics in ecological networks)

force.m is the ODEs of the plant-pollinator mutualistic model;

multivariate_normal_distribution.m is the density function of the Gaussian distribution which is used to calculate the density function of the model;

Solver.m is used to solve the ODEs, calculate the mean value and the covariance of each stable state.

minActionPath.m is used to calculate the transition paths and actions between the stable states.

calculate_sigma.m is used to calculate the covariance of one stable state.

main.m is the main function to calculate the density function of species abundance and use our dimension reduction approach to plot the landscape and transition paths of the system.

M_PL_006_tristable.xlsx is the tristable network which we analyzed in our main text.

Please run main.m to simulate our program.
