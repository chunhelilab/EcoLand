function f=force(t,x,p,r1,r2,net)
    %alphaP=p(1);alphaA=p(2);betaP=p(3);betaA=p(4);h=p(5);
    %uP=p(6);uA=p(7);kappa=p(8);betaother=p(9);r0=p(10);delta=p(11);
    f = zeros(size(net,1)+size(net,2),size(x,2));
    for ii = 1:size(net,1)
    f(ii,:) = p(1)*x(ii,:)-p(3)*x(ii,:).*x(ii,:)-p(9)*x(ii,:).*sumation(x(1:size(net,1),:))+x(ii,:).*...
                mut(x(1:size(net,1)+size(net,2),:),r1,ii)./(1+p(5).*mut(x(1:size(net,1)+size(net,2),:),r1,ii))+p(6);
    end
    
    for qq = size(net,1)+1:size(net,1)+size(net,2)
    f(qq,:) = p(2)*x(qq,:)-p(4)*x(qq,:).*x(qq,:)-p(9)*x(qq,:).*sumation(x(size(net,1)+1:size(net,1)+size(net,2),:))-p(8)*x(qq,:)...
    +x(qq,:).*mut2(x(1:size(net,1)+size(net,2),:),r2,qq)./(1+p(5)*mut2(x(1:size(net,1)+size(net,2),:),r2,qq))+p(7);
    end
end

function sumvec = sumation(X)
    sumvec = 0;
    for k = 1:size(X,1)
    sumvec = sumvec + X(k,:);
    end
end

function mutfrompolli = mut(Y,r1,ii)
    mutfrompolli = 0;
    for j = 1:size(r1,2)
    mutfrompolli = mutfrompolli+r1(ii,j).*Y(size(r1,1)+j,:);
    end
end

function mutfromplant = mut2(Z,r2,qq)
    mutfromplant = 0;
    for j = 1:size(r2,1)
    mutfromplant = mutfromplant+r2(j,qq-size(r2,1)).*Z(j,:);
    end
end
