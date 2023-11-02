function [xx,sigma,n,ycell,action,ActionVal]=Solver(cycle_index,par,d,N,net)
%%cycle_index: the number of random initial conditions to the ODEs to be solved
%%par: the parameters of the ODE
%%d: the diffusion coefficient
%%N: the dimension of the system

%%%parameters for calculating paths
params={};
params.N = 20; %The number of points in the minimum action path.
params.TMax = 50; %The time range of the minimum action path.
%Larger values can be more accurate, but can also lead to instabilities.
params.Dimension = N; %The dimension of the system.
params.c = 1e10; %The remeshing parameter, c. Larger values give greater remeshing.
params.MaxIter = 1e5; %The number of optimization steps to run between remeshing steps.
params.K= 2; %The number of total remeshing steps to do;
params.Remesh = 1; %Turn the remesh on or off. If off, the algorithm will be default perform K*MaxIter itereations.
params.q = 3; %The q parameter from the original paconstraints can be set of the path as well. This is not done in our per. Is a measure of path smoothness.
params.ub = []; %If desired, manuscript.
params.lb = []; %As above, lower bounds can be set on the path as well. This is also not done at all in our manuscript.
params.PhiInit =[]; %The default intial path between two states is a straight line. This can be modified here.
%params.ObjGrad = Parameters.JacobianPresent;%We require a Jacobian in this implementation.
params.ObjGrad=0;

xx=zeros(cycle_index,N);

gamma_0 = par(10); delta = par(11);
r1 = net*gamma_0./(sum(net')'.^delta);
r2 = net*gamma_0./(sum(net).^delta);
r1(isnan(r1)) = 0;
r2(isnan(r2)) = 0;

%%Solve odes from different initial values
for i=1:cycle_index
    x0=rand(N,1);
    [t,x]=ode45(@(t,x)force(t,x,par,r1,r2,net),[0,100],x0);
    newx=x(end,:);
    x=inf*ones(1,N);
    while norm( x(end,:)-newx(end,:) ,2 )>1e-3
        x=newx;
        [t,newx]=ode45(@(t,x)force(t,x,par,r1,r2,net),[0,1],x(end,:));
    end
    xx(i,:)=newx(end,:);
end

%%Finding the stable points
 for q=1:(cycle_index-1)
     for p=(q+1):cycle_index
         if norm(xx(q,:)-xx(p,:),'fro')<10^-1
             xx(p,:)=xx(q,:);
         end
     end
 end
stable_point=unique(xx(:,:),'rows');
n=zeros(1,2);
sigma=zeros(size(xx,1),size(xx,2)^2);
for i=1:size(stable_point,1)
    [m]=find(xx(:,2)==stable_point(i,2));
    if length(m)>=1
        disp(strcat(num2str(stable_point(i,:)),' repeat ',num2str(length(m)),' times',' the location in the row xx is' ,mat2str(m)))
    end
    n(i,1)=m(1);
    n(i,2)=length(m);
    %%%calculate the covariance of each stable state
     sig=calculate_sigma(xx(m(1),:),par,N,d,r1,r2,net)';  
     for j=1:length(m)
         sigma(m(j),:)=sig;
     end
end

%Arrange the index n from large to small by fourth elements
m=size(n,1);
if(m~=1)
    for i=1:m
        tran=n(i,:);
        flag=i;
        if( i~=m )
            for j=i+1:m
                if( xx(n(j,1),1)>xx(tran(1),1) )
                    tran=n(j,:);
                    flag=j;
                end
            end
        end
        n(flag,:)=n(i,:);
        n(i,:)=tran;
    end
end
%%stable state
SS=zeros(size(xx,2),size(n,1));
for i=1:size(n,1)
    SS(:,i)=xx(n(i,1),:)';
end

Func=@(x)force(1,x,par,r1,r2,net);  %Force
dFunc=[]; %Jacobian
feasi_n=size(n,1); 
action=zeros(feasi_n,feasi_n);
ycell=cell(feasi_n,feasi_n);
ActionVal=cell(feasi_n,feasi_n);
for i=1:feasi_n
    for j=1:feasi_n
        if i==j
            action(i,j) = Inf;
        else
            Spot1 = i;
            Spot2 = j;%initial and end point
            params.init = SS(:, [Spot1,Spot2]);
            [ActionVal{i,j}, Path] = minActionPath(params,params.init,Func,dFunc);
            action(i,j) = ActionVal{i,j}(end);
            ycell(i,j) = {[SS(:,i) Path SS(:,j)]};
        end
    end
end
end
