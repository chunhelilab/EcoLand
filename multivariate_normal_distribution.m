function z=multivariate_normal_distribution(x,x0,sigma,n)
z=1/((2*pi)^(n/2)*det(sigma)^(1/2))*exp(1)^(-0.5*(x-x0)'*sigma^(-1)*(x-x0));