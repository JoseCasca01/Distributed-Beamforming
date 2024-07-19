close all;
clear;
clc;

f = 6e9; % 6 GHz
c = 3e8;
lambda = c/f; 
flag = true;

while flag
    N = input('Number of nodes: ');
    if (mod(N,1)==0) && (N > 0)
        flag = false; 
    end
end

%Define a rectangular field
fieldx = 5; %0 to fieldx
fieldy = 5; %0 to fieldy

groupx = 2*lambda;
groupy = 2*lambda;

%Base station coordinates
BS = [fieldx/2,fieldy/2,2*(10*lambda)^2/lambda]; %FarField
%BS = [fieldx/2,fieldy/2,20*lambda];
if BS(3) < 1
    BS(3) = 1;
elseif BS(3)>10
    BS(3) = 10;
end
%Define the number of sensors per group
%5 groups
sensorsPos = zeros(N*5,4); %rows - sensors; collums - coordenates x and y

%The groups need to have sensors close to each other
%Groups dimension 1.5 by 1.5
%Define groups limits 
%sensorsPos -> 0.5 to 2 on X and 0.5 to 2 on Y

%sensorsPos(1:end,1) = 0.5 + 1.5 * rand(size(sensorsPos,1),1);
%sensorsPos(1:end,2) = 0.5 + 1.5 * rand(size(sensorsPos,1),1);

deslocamentox = [0 1.2 1.4 4 3];
deslocamentoy = [1.2 4 0.5 1.4 3];

%Random sensors in all the field
for i = 1:length(deslocamentox)
    for j = 1:N
        sensorsPos(N*(i-1)+j,1) = deslocamentox(i) + groupx * rand(1,1);
        sensorsPos(N*(i-1)+j,2) = deslocamentoy(i) + groupy * rand(1,1);
        sensorsPos(N*(i-1)+j,4) = i;
    end
end

figure(1);
plot3(BS(1),BS(2),BS(3),'rX'), hold on;
for i = 1:length(deslocamentox)
    plot3(sensorsPos(sensorsPos(:,4)==i,1),sensorsPos(sensorsPos(:,4)==i,2),sensorsPos(sensorsPos(:,4)==i,3),'O');
end
title(['Field with ', num2str(N), ' nodes']);
ylabel('yfield (m)');
xlabel('xfield (m)');
zlabel('height (m)')
legend('Base Station', 'Group 1', 'Group 2', 'Group 3', 'Group 4', 'Group 5');
axis([0, fieldy, 0, fieldx]);

%% Fases ajustada para cada sensor

R = zeros(5,N);

for i = 1:length(deslocamentox)
    R(i,:) = distance(sensorsPos(sensorsPos(:,4)==i,:),BS);
end
t = linspace(0,1.2e-7,1000000);
totalsignal = zeros(1,length(t));

[t, adjustsignals, adjustsignals_Loss, receivedsignal,receivedsignal_Loss, PropagationLoss] = receptorgroups(sensorsPos,R,max(max(R)),f,c,N,lambda,1);

figure;
plot(t,adjustsignals(:,:)),hold on;
title('Received Signals');
xlabel('t(s)');
ylabel('Amplitude');
axis([3.32e-8 3.50e-8 -1 1]);

figure; 
plot(t,receivedsignal_Loss),hold on;
title('Received Signal');
xlabel('t(s)');
ylabel('Amplitude');
axis([3.32e-8 3.50e-8 min(receivedsignal_Loss) max(receivedsignal_Loss)]);

maxamp_exactphase = max(receivedsignal_Loss);

%% Fase ajustada para cada grupo


[t, adjustsignals, adjustsignals_Loss, receivedsignal,receivedsignal_Loss, PropagationLoss] = receptorgroupapprox(sensorsPos,BS,f,c,N,lambda,1);

%%
cor = ['b' 'r' 'k' 'g' 'y'];
figure;
for i = 1:5
    plot(t,adjustsignals(1+20*(i-1):20+20*(i-1),:),cor(i)),hold on;
end
%title('Received Signals');
xlabel('t(s)');
ylabel('Amplitude');
axis([3.46e-8 3.50e-8 -1 1]);

figure;
plot(t,receivedsignal_Loss),hold on;
%title('Received Signal');
xlabel('t(s)');
ylabel('Amplitude');
axis([3.32e-8 3.50e-8 min(receivedsignal_Loss) max(receivedsignal_Loss)]);

maxamp_approxphase = max(receivedsignal_Loss);