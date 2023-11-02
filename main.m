clear
cycle_index = 300;  %% The number of random initial conditions to the ODEs to be solved
alphaP = 0.3; alphaA = 0.3; betaP = 1; betaA = 1; h = 0.2; uP = 1e-3; uA = 1e-3; 
kappa = 1.1; betaother = 0.01; gamma_0 = 1; delta = 0.5;
par = [alphaP,alphaA,betaP,betaA,h,uP,uA,kappa,betaother,gamma_0,delta]; %The parameters of the plant-pollinator network
d = 0.0005;  %% The diffusion coefficient    

[net] = xlsread('M_PL_006_tristable.xlsx'); %tristable system

N=size(net,1)+size(net,2); %% The dimension of the system

    %% Solve the ODEs, calculate the paths and actions;
[xx,sigma,n,ycell,action,ActionVal]=Solver(cycle_index,par,d,N,net);

    index=size(n,1);  %% The number of the stable states
    alpha=zeros(index,1);  %% The weight of the stable states
    sigma0=cell(index,1);  %% The covariance of the Gaussian density function
    mu=zeros(index,N);  %% The mean value of the Gaussian density function

    for i=1:index
        %The mean value of each stable state
        mu(i,:)=xx(n(i,1),:); 
        %The covariance of each stable state
        sigma0{i}=reshape(sigma(n(i,1),:),N,N)';  
        %The weight of each stable state
        alpha(i)=n(i,2)/sum(n(:,2)); 
    end
      %% DRL
    %Calculate the mean value of plants and pollinators
    Mu1=0;
    for i=1:index
        Mu1=Mu1+alpha(i)*mu(i,1:size(net,1));
    end
    
    Mu2=0;
    for i=1:index
        Mu2=Mu2+alpha(i)*mu(i,size(net,1)+1:end);
    end
    
    %Calculate the covariance of plants and pollinators
    Sigma1=-Mu1'*Mu1; 
    for i=1:index
        Sigma1=Sigma1+alpha(i)*(sigma0{i}(1:size(net,1),1:size(net,1))+mu(i,1:size(net,1))'*mu(i,1:size(net,1)));
    end
    
    Sigma2=-Mu2'*Mu2; 
    for i=1:index
        Sigma2=Sigma2+alpha(i)*(sigma0{i}(size(net,1)+1:end,size(net,1)+1:end)+mu(i,size(net,1)+1:end)'*mu(i,size(net,1)+1:end));
    end

    %Calculate the eigenvalues and eigenvectors of the covariance
    [V1,D1] = eigs(Sigma1,1);
    [V2,D2] = eigs(Sigma2,1);
    V = [V1;zeros(size(V2,1),1);zeros(size(V1,1),1);V2];
    V = reshape(V,size(net,1)+size(net,2),2);
    
    if sign(V(:,1)'*ones(N,1))<0
        V(:,1)=-V(:,1);
    end
    if sign(V(:,2)'*ones(N,1))<0
        V(:,2)=-V(:,2);
    end

    %%%Calculate the covariance and mean value after dimension reduction
    sigma0_pca=cell(index,1);
    mu_pca=zeros(index,2);
    for i=1:index
        mu_pca(i,:)=V'*mu(i,:)';
        sigma0_pca{i}=V'*sigma0{i}*V;
    end

    %% plot the landscape
    y_max=[6,3]; %% Range of the landscape
    y_min=[-1,-1];
    step=(y_max-y_min)/100; %% Length of the step
    [a1,a2]=meshgrid(y_min(1):step(1):y_max(1),y_min(2):step(2):y_max(2)); %% Grid
    [s1,s2]=size(a1);
    P=zeros(s1,s2);
    z=zeros(s1,s2);
    for kk=1:index
        sig=sigma0_pca{kk};
        x_wen=mu_pca(kk,:);
        for i=1:s1
            for j=1:s2
                z(i,j)=multivariate_normal_distribution([a1(i,j);a2(i,j)],x_wen',sig,2);  %% Normal distribution
            end
        end

        P=P+z*alpha(kk);
    end
    P=real(P);
    P=P/sum(sum(P));
    figure(time);
    surf(a1,a2,-log(max(P,10^-100)));   %% Plot landscape
    shading interp
    xlabel('PC1')
    ylabel('PC2')
    zlabel('U')
    axis([-1 6 -1 3 0 240])

    for i=1:size(n,1)
        A(i)=floor((mu_pca(i,1)-y_min(1))/step(1))+1;
        B(i)=floor((mu_pca(i,2)-y_min(2))/step(2))+1;
    end
    hold on

    %Plot the grid
    for i=1:floor(size(a1,1)/4)
        plot3(a1(4*i-1,:),a2(4*i-1,:),-log(max(P(4*i-1,:),10^-100)),'Color',[0.4 0.4 0.4],'LineWidth',0.01);
    end
    for i=1:floor(size(a1,2)/4)
        plot3(a1(:,4*i-1),a2(:,4*i-1),-log(max(P(:,4*i-1),10^-100)),'Color',[0.4 0.4 0.4],'LineWidth',0.01);
    end

    %% Plot the paths
    %Calculate the paths after dimension reduction
        y12=V'*ycell{1,2};y13=V'*ycell{1,3};
        y21=V'*ycell{2,1};y23=V'*ycell{2,3};
        y31=V'*ycell{3,1};y32=V'*ycell{3,2};
        view([-25,75])
        hold on
        z3path=griddata(a1,a2,-log(max(P,10^-100)),y12(1,:),y12(2,:));
        plot3(y12(1,:),y12(2,:),z3path+3,'w','LineWidth',2);
        z3path=griddata(a1,a2,-log(max(P,10^-100)),y21(1,:),y21(2,:));
        plot3(y21(1,:),y21(2,:),z3path+3,'Color',[0.85,0.43,0.83],'LineWidth',2);
        z3path=griddata(a1,a2,-log(max(P,10^-100)),y13(1,:),y13(2,:));
        plot3(y13(1,:),y13(2,:),z3path+3,'w','LineWidth',2);
        z3path=griddata(a1,a2,-log(max(P,10^-100)),y31(1,:),y31(2,:));
        plot3(y31(1,:),y31(2,:),z3path+3,'Color',[0.85,0.43,0.83],'LineWidth',2);
        z3path=griddata(a1,a2,-log(max(P,10^-100)),y23(1,:),y23(2,:));
        plot3(y23(1,:),y23(2,:),z3path+3,'w','LineWidth',2);
       z3path=griddata(a1,a2,-log(max(P,10^-100)),y32(1,:),y32(2,:));
        plot3(y32(1,:),y32(2,:),z3path+3,'Color',[0.85,0.43,0.83],'LineWidth',2);

        view([-25 75])
        set(gcf,'outerposition', [100 100 800 650]);