clear all
clc
PosTarget = [-75,75]; %true position of the target
[Ntx,Nrx,Mtx,Nue,gammac,nvariance,lambdac,Dsn,PosAPs,Posue,PosRx,distTx,distRx,Txbeamsteering,Rxbeamsteering,alphark,sigmarcs,H,sigma2,W]=ParametersCont(PosTarget);
xcont = -475:50:475; %center of the search area cells
ycont = 475:-50:-475;
%n = length(xcont) ;
%m = length(ycont) ; 
[Xcont,Ycont] = meshgrid(xcont,ycont) ; 
x =-500:50:500; %Entire search area
y = 500:-50:-500;
n = length(x) ;
m = length(y) ; 
[X,Y] = meshgrid(x,y) ; 
pd = 0;
pa = 0;
Poscont = zeros(100,2);
Mue = [2 4 8 16]; %modify this part indication the number of antennas at the UE or modify accordingly to evaluate other parameter.
tic
%for itery = 1:length(Mue)
itery = 3;
parfor MC = 1:100 %obtaining the covariance matrix with the EM algorithm
    Rall(:,:,:,:,MC)= Em(Dsn,nvariance,Ntx,Nue,Mue(itery),H((1+(itery-1)*(Mue(itery))):itery*(Mue(itery)),:,:,MC),sigma2,W(:,(1+(itery-1)*(Nue+1)):itery*(Nue+1),:,MC),Mtx);
end
    %Rtotal(:,:,:,:,(1+(itery-1)*100):itery*100) = Rall;
%end
%figure
parfor MC = 1:100 %verification for each MC iteration of the probability of detection
    angmax = zeros(1,Ntx);
if Ntx ==8 %interval of the steering vectors 
angsAPs = [0,pi/2;-pi,-pi/2;-pi/2,0;pi/2,pi;-pi/2,pi/2;-pi/2,pi/2;pi/2,-pi/2;-pi/2,pi/2];
elseif Ntx == 4
angsAPs = [-pi/2,pi/2;-pi/2,pi/2;-pi/2,pi/2;-pi/2,pi/2];
else 
angsAPs = [pi/2,-pi/2;pi/2,-pi/2];
end
ptotal = zeros(901,128);
pt = zeros(1,901,Ntx);
for k =1:Ntx
    %if abs(PosAPs(k,1)) == abs(PosAPs(k,2))
        thetal = linspace(angsAPs(k,1),angsAPs(k,2),901); %angle grid sample
        a = (1i*2*pi*(1/2)*(1:1:(Mtx-1)).*sin(thetal'))'; %beam steering
        %atargets = (i*2*pi*(1/2)*(1:1:(M-1)).*sin([theta2 theta1 theta3]'))';
        asteer = [ones(1,length(thetal));exp(a)];
        py = zeros(1,length(thetal));
        for temp = 1:128
        for i=1:length(thetal)
            py(i) = real(asteer(:,i)'*Rall(:,:,temp,k,MC)*asteer(:,i)); %beampattern
        end
        ptotal(:,temp) = py;    
        end
        pt(:,:,k) = sum(ptotal,2)'; %beampattern for all the APs
end
contgrid = zeros(length(x)-1,length(y)-1,Ntx);
[maxim,ind] = max(pt); %get the index and value of the maximum point of the beampattterns
xreta = -500:0.001:500; 
yreta = zeros(1,length(xreta),Ntx);
angcont = [0,pi/2;-pi,-pi/2;-pi/2,0;pi/2,pi;0,pi;0,-pi;pi/2,-pi/2;-pi/2,pi/2];
for k = 6:Ntx-2
    thetal = linspace(angcont(k,1),angcont(k,2),901); %angle grid sample
    angmax(k) = thetal(ind(:,:,k));
    yreta(:,:,k) = tan(angmax(k)).*(xreta-PosAPs(k,1))+PosAPs(k,2); %calculate the direction line
for i = 1:m-1
  for j = 1:n-1
    % Get the points of each grid 
    P = [X(i,j) Y(i,j) ; 
         X(i,j+1) Y(i,j+1) ; 
         X(i+1,j+1) Y(i+1,j+1) ;
         X(i+1,j) Y(i+1,j)] ; 
         
    idx = any(inpolygon(xreta,yreta(:,:,k),P(:,1),P(:,2)));
    contgrid(i,j,k) = idx;
    %iwant = [x(idx) y(idx) z(idx)] ; 
    %plot(x(idx),y(idx),'.')
    %drawnow
   end
end
end

contTotal = (sum(contgrid,3)); %obtain a matrix with the count of cells that are crossed by the lines
maximum = max(max(contTotal)); %verify the line with the maximum number of lines crossed.
[xc,yc]=find(contTotal==maximum);
if length(xc)>1
    aux = randi(length(xc));
    xf = xc(aux);
    yf = yc(aux);
    Poscont(MC,:) = [Xcont(xf,yf),Ycont(xf,yf)];
else
    xf = xc;
    yf = yc;
    Poscont(MC,:) = [Xcont(xf,yf),Ycont(xf,yf)];
end
if Poscont(MC,:) ==-PosTarget %calculate the probability of detection
    pd = pd + 1;
end
end
%Posconttotal = [Posconttotal Poscont];
toc
save Posconts1knownAPsm81.mat Poscont