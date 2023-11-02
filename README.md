# EcoLand
EcoLand is a pipeline to simulate the stochastic dynamics of the ecological systems (Jinchao Lv, Jin Wang* and Chunhe Li*, Landscape quantifies the intermediate state and transition dynamics in ecological networks)

force.m is the ODEs of the plant-pollinator mutualistic model;

multivariate_normal_distribution is the density function of the Gaussian distribution which is used to calculate the density function of the model;

Solver.m is used to solve the ODEs, calculate the mean value and the covariance of each stable state and calculate the transition paths and actions between the stable states.

calculate_sigma.m is used to calculate the covariance of one stable state.

main.m is the main function to calculate the density function of species abundance and use our dimension reduction approach to plot the landscape and transition paths of the system.

Please run main.m to simulate our program.
