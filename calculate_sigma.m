function [sig]=calculate_sigma(xx,par,kk,d,r1,r2,net)
alphaP=par(1);alphaA=par(2);betaP=par(3);betaA=par(4);h=par(5);uP=par(6);
uA=par(7);kappa=par(8);betaother=par(9);r0=par(10);epsilon=par(11);

syms P
for i=1:size(net,1)
    P(i)=sym(['P' num2str(i)]);
end

syms A
for i=1:size(net,2)
    A(i)=sym(['A' num2str(i)]);
end

for i=1:size(net,1)
    mm(i)=xx(i);
end

for i=1:size(net,2)
    nn(i)=xx(i+size(net,1));
end
   
Ajac=jacobian([alphaP*P'-betaP*P'.*P'-betaother*P'*sum(P)+P'.*r1*A'./(1+h*r1*A')+uP;...
     alphaA*A'-betaA*A'.*A'-betaother*A'*sum(A)-kappa*A'+A'.*r2'*P'./(1+h*r2'*P')+uA],...
     [P,A]);

Ajac=double(subs(Ajac,[P,A],[mm,nn]));

% A*sigma+sigma*A'+2D

P=zeros(kk^2,kk^2);  %coefficient matrix

%%the initial of coeffiicient matrix
for i=0:(kk-1)
    P(i*kk+1:i*kk+kk,i*kk+1:i*kk+kk)=P(i*kk+1:i*kk+kk,i*kk+1:i*kk+kk)+Ajac;
end

for m=0:kk-1
    for i=1:kk
        for j=1:kk
            P(m*kk+i,(j-1)*kk+i)=P(m*kk+i,(j-1)*kk+i)+Ajac(m+1,j);
        end
    end
end

B=zeros(kk^2,1);
for i=1:kk
    B((i-1)*kk+i)=-2*d;
end


sig=P\B;

end