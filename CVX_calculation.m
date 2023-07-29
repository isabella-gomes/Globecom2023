function [Wtotal,obj]=CVX_calculation(PosTarget,Ntx,Nrx,Mtx,Nue,Mue,H,gammac,nvariance)

[lambdac,Dsn,PosAPs,Posue,PosRx,distUe,distTx,distRx,Txbeamsteering,Rxbeamsteering,alphark,sigmarcs,Hno]=Parameters(Ntx,Nrx,Mtx,Nue,Mue,PosTarget);
sigma2 = abs(nvariance.*randn(1,1));
obj = -20 + (20 + 20)*rand(Mue,Nue); %initial value for the receive beamforming at the users
W1 = -20 + (20 + 20)*rand(Mtx,Nue+1,Ntx); %initial value for the transmit beamformer with index K for the optimization of K==J
W2 = W1; %initial value for the transmit beamformer with index K for the optimization of K~=J
%wjtil = rand(Mtx,Nue+1,Ntx); %initial value for the transmit beamformer with index J
%obj = zeros(Mue,Nue);                                                      
Txbeam = permute(Txbeamsteering,[3 2 1]);
Rxbeam = permute(Rxbeamsteering,[3 2 1]);
covalpha = (conj((sigmarcs).*(randn(1,1)+1i*randn(1,1)))*((sigmarcs).*(randn(1,1)+1i*randn(1,1))));
Ptotal = 10^(50/10)/1000;
convergence = 0;
Miter = 0;
%A = zeros(Nue+1,Nue+1); 
while convergence==0
    rbf = obj;
    Wtil1 = W1+W2;
    Wtil2 = W1+W2;
    %A1 = zeros(Nue+1,Nue+1);
    %A2 = zeros(Nue+1,Nue+1);
    cvx_clear;    
    cvx_begin quiet%calculating the case where K == J 
        variable W1(Mtx,Nue+1,Ntx) 
        variable a1(Ntx,Nue) nonnegative
        variable b1(Ntx,Nue) nonnegative
        expression A1(Nue+1,Nue+1);
        expression Atil1(1,1);
        for r = 1:Nrx
            for k = 1:Ntx
                Asteer = conj(Txbeam(:,k))*Txbeam(:,k)';
                constante = (sqrt((power(lambdac,2)/(power(4.*pi,3).*power(distTx(k),2).*power(distRx(r),2))).*(power(lambdac,2)./(power(4.*pi,3).*power(distTx(k),2).*power(distRx(r),2))))*Rxbeam(:,r,:)'*covalpha*Rxbeam(:,r,:));
                A1 = A1 + real(constante.*(Wtil1(:,:,k)'*Asteer*Wtil1(:,:,k)+2*Wtil1(:,:,k)'*Asteer*(W1(:,:,k)-Wtil1(:,:,k))));
            end
        end
        for n = 1:length(Dsn)
            Atil1 = Atil1 + Dsn(:,n)'*A1*Dsn(:,n);
        end
        Af1 = real(Atil1)/(Mtx*Nrx*sigma2);
        %obj = real(obj * sigmasq_sens);

        maximize Af1
         subject to
            for k = 1:Ntx
                den = 0;
                 for i = 1:Nue
                    real((power(rbf(:,i)'*H(:,(1+(k-1)*Mtx):k*Mtx,i)*Wtil1(:,i,k),2) + 2*Wtil1(:,i,k)'*H(:,(1+(k-1)*Mtx):k*Mtx,i)'*rbf(:,i)*rbf(:,i)'*H(:,(1+(k-1)*Mtx):k*Mtx,i)*(W1(:,i,k)-Wtil1(:,i,k))))>=a1(k,i);
                    Waux = W1(:,:,k);
                    Waux(:,i,:)=[];
                    b1(k,i)>=(sum(power(real(rbf(:,i)'*H(:,(1+(k-1)*Mtx):k*Mtx,i)*Waux),2))+sigma2.*norm(rbf(:,i),2));
                    a1>=gammac*b1(k,i);
                 end
                norm(W1(:,:,k))<=Ptotal;
            end
    cvx_end
    C1 = cvx_optval;

    cvx_clear;    
    cvx_begin quiet%calculating the case where K ~= J 
        variable W2(Mtx,Nue+1,Ntx) 
        variable a2(Ntx,Nue) nonnegative
        variable b2(Ntx,Nue) nonnegative
        expression A2(Nue+1,Nue+1);
        expression Atil2(1,1);
        for r = 1:Nrx
            for k = 1:Ntx
                Waux = Wtil2;
                Waux(:,:,k)=[];
                Txaux = Txbeam;
                Txaux(:,k)=[];
                for j = 1:Ntx-1
                    A2 = A2 + real(sqrt((power(lambdac,2)/(power(4.*pi,3).*power(distTx(k),2).*power(distRx(r),2))).*(power(lambdac,2)./(power(4.*pi,3).*power(distTx(j),2).*power(distRx(r),2))))*Rxbeam(:,r,:)'*covalpha*Rxbeam(:,r,:)*W2(:,:,k)'*conj(Txbeam(:,k))*Txaux(:,j)'*Waux(:,:,j));
                end
            end
        end
        for n = 1:length(Dsn)
            Atil2 = Atil2 + Dsn(:,n)'*A2*Dsn(:,n);
        end
        Af2 = real(Atil2)/(Mtx*Nrx*sigma2);
        %obj = real(obj * sigmasq_sens);

        maximize Af2
         subject to
            for k = 1:Ntx
                den = 0;
                 for i = 1:Nue
                    real((power(rbf(:,i)'*H(:,(1+(k-1)*Mtx):k*Mtx,i)*Wtil2(:,i,k),2) + 2*Wtil2(:,i,k)'*H(:,(1+(k-1)*Mtx):k*Mtx,i)'*rbf(:,i)*rbf(:,i)'*H(:,(1+(k-1)*Mtx):k*Mtx,i)*(W2(:,i,k)-Wtil2(:,i,k))))>=a2(k,i);
                    Waux = W2(:,:,k);
                    Waux(:,i,:)=[];
                    b2(k,i)>=(sum(power(real(rbf(:,i)'*H(:,(1+(k-1)*Mtx):k*Mtx,i)*Waux),2))+sigma2.*norm(rbf(:,i),2));
                    a2(k,i)>=gammac*b2(k,i);
                 end
                norm(W2(:,:,k))<=Ptotal;
            end
    cvx_end
    C2 = cvx_optval;
    Wtotal = W1+W2;
    Waux = Wtil1+Wtil2;
    for k = 1:Ntx
    if norm(Wtotal(:,:,k)-Waux(:,:,k),'fro')/norm(Wtotal(:,:,k),'fro')<= 0.5
            convergence=1;
            Wotimo=Wtotal;
    end
    end
     sumWc = zeros(Mue, Mue, Nue);
     wi = zeros(Mtx*Ntx,Nue+1);
    for k = 1:Ntx
        for i = 1:Nue 
            Wa = Wtotal(:,:,k);
            Wa(:,i,:)=[];
            for j = 1:Nue-1
                sumWc(:,:,i) = sumWc(:,:,i) + real(H(:,(1+(k-1)*Mtx):k*Mtx,i)*Wa(:,j,:)*Wa(:,j,:)'*H(:,(1+(k-1)*Mtx):k*Mtx,i)');
            end
        end
        wi((1+(k-1)*Mtx):k*Mtx,:) = Wtotal(:,:,k);
    end

    for i = 1:Nue
         obj(:,i) = real(inv(sumWc(:,:,i)+sigma2*eye(Mue))*H(:,:,i)*wi(:,i));
    end 
    %display(beam)
    normu = norm(obj-rbf)/norm(obj);
    if isnan(normu)
        obj = -20 + (20 + 20)*rand(Mue,Nue); %initial value for the receive beamforming at the users
        W1 = -20 + (20 + 20)*rand(Mtx,Nue+1,Ntx); %initial value for the transmit beamformer with index K for the optimization of K==J
        W2 = W1; %initial value for the transmit beamformer with index K for the optimization of K~=J
        Miter=Miter-1;
    end
    if Miter >=10
        if C1 == inf && C2 == inf
            continue
        else
            break
        end
    end
    Miter = Miter+1;
    if norm(Wtotal -Wtil1)/norm(Wtotal)<0.1 || normu<0.1
        convergence = 1;
    end
    
end
end