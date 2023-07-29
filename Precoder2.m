clear all
tic
xtarget = -75;
ytarget = -xtarget;
%PosTarget = [xtarget',-125*ones(length(xtarget),1)];
Ntx = 8; %number of transmit APs
Nrx = 4; %number of receive APs
Mtx = 64; %number of antennas at all the APs
Nue = 4; %number of users
Mue = [2 4 8 16]; %number of antennas at the users
MC = 100;
%H4users = load('HUe.mat');
%H = H4users.H;
gammac = 10^(3/10); %threshold for the communication
nvariance = 10^(-94/100)/1000;
%[lambdac,Dsn,PosAPs,Posue,PosRx,distUe,distTx,distRx,Txbeamsteering,Rxbeamsteering,alphark,sigmarcs,H]=Parameters(Ntx,Nrx,Mtx,Nue,Mue,[xtarget(1,1),xtarget(1,1)]);
Wx = zeros(Mtx,Nue+1,Ntx,MC);
%W = zeros(Mtx,(Nue+1)*length(xtarget),Ntx,MC);

for itery=1:length(Mue)
ux = zeros(Mue(itery),(Nue),MC);
HMCx = zeros(Mue(itery),Mtx*Ntx,Nue,MC);
parfor iterx=1:MC
[lambdac,Dsn,PosAPs,Posue,PosRx,distUe,distTx,distRx,Txbeamsteering,Rxbeamsteering,alphark,sigmarcs,H]=Parameters(Ntx,Nrx,Mtx,Nue,Mue(itery),[xtarget,ytarget]);
[Wtotal,obj]=CVX_calculation([xtarget,ytarget],Ntx,Nrx,Mtx,Nue,Mue(itery),H,gammac,nvariance);
Wx(:,:,:,iterx) = Wtotal;
ux(:,:,iterx) = obj;
HMCx(:,:,:,iterx) = H;
%iterx
end
HMC((1+(itery-1)*(Mue(itery))):itery*(Mue(itery)),:,:,:)=HMCx;
W(:,(1+(itery-1)*(Nue+1)):itery*(Nue+1),:,:)=Wx;
u((1+(itery-1)*(Mue(itery))):itery*(Mue(itery)),:,:)=ux;
itery
end
toc
save HMC2.mat HMC
save WMC2.mat W
save uMC2.mat u
