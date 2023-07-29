function [lambdac,Dsn,PosAPs,Posue,PosRx,distUe,distTx,distRx,Txbeamsteering,Rxbeamsteering,alphark,sigmarcs,H]=Parameters(Ntx,Nrx,Mtx,Nue,Mue,PosTarget)
lambdac = 1.9e9; %carrier frequency
[cn,sn]=signalgeneration(Nue,1); %generating symbols for radar and communication
Dsn = [sn;cn]; %diagonal matrix containing the sensing and communication symbols
%PosAPs = -500 + (500 + 500).*rand(Ntx,2); %Position of transmit APs
%Posue = -500 + (500 + 500).*rand(Nue,2); %Position of Users
PosAPs = [-500,-500;500,500;-500,500;500,-500;0,-500;0,500;-500,0;500,0]; %location of APs considering 8 Ntx
%PosAPs = [0,-500;0,500;-500,0;500,0]; %location of APs considering 4 Ntx
%PosAPs = [0,-500;0,500]; %location of APs considering 2 Ntx
Posue = [300,300;-300,-300;-300,300;300,-300];
%PosTarget = [25,-350];
%PosRx = -15 + (15 + 15).*rand(Nrx,2);
PosRx = [250,250;-250,250;250,-250;-250,-250]; %location of Receive APs for Nrx= 4
%PosRx = [-250,250;250,-250]; %location of Receive APs for Nrx = 2
%PosRxaux = repmat(PosRx,Ntx-Nrx,1);
%Posueaux = repmat(Posue,Ntx,1);
%PosAPaux = repmat(PosAPs,Nue,1);
%Rrk = corr2(PosRxaux,PosAPs);
%distTR = sqrt((PosRxaux(:,1)-PosAPs(:,1)).^2+(PosRxaux(:,2)-PosAPs(:,2)).^2);

%distUe =  sqrt((Posueaux(:,1)-PosAPs(:,1)).^2+(Posueaux(:,2)-PosAPs(:,2)).^2);
for i=1:Ntx
    for j = 1:Nue
        distUe(i,j) =  sqrt((Posue(j,1)-PosAPs(i,1)).^2+(Posue(j,2)-PosAPs(i,2)).^2);
    end
end

[angleTx, distTx] = compute_angle_dist(PosAPs, PosTarget); %azimuth angle calculation with transmit APs
[angleRx, distRx] = compute_angle_dist(PosRx, PosTarget); %azimuth angle calculation with receive APs
Txbeamsteering = beamsteering(angleTx.', Mtx); % (Ntarget x Ntx x Mtx) %transmit beamsteering
Rxbeamsteering = beamsteering(angleRx.', Mtx); % (Ntarget x Nrx x Mtx) %receive beamsteering
sigmarcs = 10^(10/10);
alphark = (sigmarcs).*(randn(1,1)+1i*randn(1,1));
%betark = lambdac^2/((4*pi)^3*

% %Channel coefficient users 
phi = 3;
for i = 1:Ntx
    Haux = (sqrt(distUe(i,:).^(-phi)/2))'.*(randn(Nue,Mtx,Mue)+1i*randn(Nue,Mtx,Mue)); %Rayleigh modelling
    Husers = permute(Haux,[3 2 1]);
    H(:,(1+(i-1)*Mtx):i*Mtx,:) = Husers; %Number of receive antennas X Number of transmit antennas x Number of BS X Number of users
end
end