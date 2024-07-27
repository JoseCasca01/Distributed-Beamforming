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
fieldx = 10*lambda; %0 to fieldx
fieldy = 10*lambda; %0 to fieldy

%Base station coordinates
BS = [fieldx/2,fieldy/2,2*fieldx^2/lambda]; %FarField
%BS = [fieldx/2,fieldy/2,20*lambda];
if BS(3) < 1
    BS(3) = 1;
elseif BS(3)>10
    BS(3) = 10;
end
%Define the number of sensors per group
sensorsPos = zeros(N,3); %rows - sensors; collums - coordenates x and y

%The groups need to have sensors close to each other
%Groups dimension 1.5 by 1.5
%Define groups limits 
%sensorsPos -> 0.5 to 2 on X and 0.5 to 2 on Y

%sensorsPos(1:end,1) = 0.5 + 1.5 * rand(size(sensorsPos,1),1);
%sensorsPos(1:end,2) = 0.5 + 1.5 * rand(size(sensorsPos,1),1);

%Random sensors in all the field
sensorsPos(1:end,1) = fieldx * rand(size(sensorsPos,1),1);
sensorsPos(1:end,2) = fieldy * rand(size(sensorsPos,1),1);

figure(1);
plot3(BS(1),BS(2),BS(3),'rX'), hold on;
plot3(sensorsPos(:,1),sensorsPos(:,2),sensorsPos(:,3),'bO');
title(['Field with ', num2str(N), ' nodes']);
ylabel('yfield (m)');
xlabel('xfield (m)');
legend('Base Station', 'Nodes');
axis([0, fieldy, 0, fieldx]);

%% Absolute position known
R = distance(sensorsPos,BS);
val = receptor(R,f,c,N,lambda,1);
%% Gaussian Error in sensors' position

figures = 1; %0 do not plot graphs / 1 plot graphs
R = distance(sensorsPos,BS);
if figures ~= 0
    receptor_SensorposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,0.2*lambda,figures);
else
    val = zeros(1,1000);
    for i = 1:length(val)
        val(i)=receptor_SensorposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,0.2*lambda,figures);
    end
    disp(mean(val));
end

%% Gaussian Error in Base Station's (BS) position

figures = 1; %0 do not plot graphs / 1 plot graphs
R = distance(sensorsPos,BS);

if figures ~= 0
    receptor_BSposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,0.2*lambda,figures);
else
    val = zeros(1,1000);
    for i = 1:length(val)
        val(i) = receptor_BSposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,0.2*lambda,figures);
    end
    disp(mean(val));
end

%% Ideal Signal Amplitude vs Sensors' position Erro
z=1;

R = distance(sensorsPos,BS);

val = receptor(R,f,c,N,lambda,0);

rounds = 5;
variance = 0.025:0.025:0.5;
Median_y =  zeros(1,length(variance));
Median_x = zeros(1,length(variance));
valtests = zeros(rounds,length(variance));

for m = 1:1
    for i = 1:length(variance)
        for j = 1:rounds*m
            valtests(j,i) = receptor_SensorposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,variance(i)*lambda,0);
        end
    end


    valnormalized = valtests/val;

    figure;
    hd = boxplot(valnormalized,variance);
    hold on;
    xlabel('Variance (wavelength)');
    ylabel('Normalized received signal amplitude');
    for i=6:7:7*length(variance)-1
        Median_x(z) = sum(get(hd(i),'XData'))/2;
        Median_y(z) = sum(get(hd(i),'YData'))/2;
        z=z+1;
    end
    plot(Median_x,Median_y,'LineWidth',3,'Color',[1.00,0.41,0.16]);
    legend('Average');
end

%% Ideal Signal Amplitude vs BS position Error
clear variance
clear valnormalized

z=1;

R = distance(sensorsPos,BS);

val = receptor(R,f,c,N,lambda,0);

rounds = 500;

%variance1(1,:) = 0.025:0.05:3;
%variance2(1,:) = 3.025:0.05:6;
%variance3(1,:) = 0.025:0.05:0.75;
%variance4(1,:) = 2:0.2:3;
variance = 3.025:0.05:4;
Median_y =  zeros(1,length(variance));
Median_x = zeros(1,length(variance));
valtests = zeros(rounds,length(variance));
  
for m = 1:1
    for i = 1:length(variance)
        for j = 1:rounds
            [valtests(j,i),~] = receptor_BSposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,variance(i)*lambda,0);
        end
    end


    valnormalized = valtests/val;

    figure;
    hd = boxplot(valnormalized,variance);
    hold on;
    xlabel('Variance (wavelength)');
    ylabel('Normalized received signal amplitude');
    for i=6:7:7*length(variance)-1
        Median_x(z) = sum(get(hd(i),'XData'))/2;
        Median_y(z) = sum(get(hd(i),'YData'))/2;
        z=z+1;
    end
    plot(Median_x,Median_y,'LineWidth',3,'Color',[1.00,0.41,0.16]);
    legend('Average');
end

%%
h=gco;
set(h,'XTick',1:4:60);
set(h,'XTickLabel',variance(1:4:60));
grid on;
% 0.835069444444445,0.898115397425656,0.069791666666667,0.027090694935218
%% Ideal Array Factor vs Sensors' position Error
clear variance
clear valnormalized
    
z = 1;

R = distance(sensorsPos,BS);

AFIdeal = N; %Array Factor (AF) for correct positions

rounds = 300;
variance = 0.05:0.025:0.5;
Median_y =  zeros(1,length(variance));
Median_x = zeros(1,length(variance));

for m = 1:1
    AFtest = zeros(rounds,length(variance));
    for i = 1:length(variance)
        for j = 1:rounds*m
            [~, Rerror] = receptor_SensorposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,variance(i)*lambda,0);
            AFtest(j,i) = sum(exp(1j*2*pi/lambda*(R-Rerror)));
        end
    end
    AFnormalized = abs(AFtest/AFIdeal);

    figure;
    hd = boxplot(AFnormalized,variance);
    hold on;
    xlabel('Variance (wavelength)');
    ylabel('Normalized Array Factor');
    for i=6:7:7*length(variance)-1
        Median_x(z) = sum(get(hd(i),'XData'))/2;
        Median_y(z) = sum(get(hd(i),'YData'))/2;
        z=z+1;
    end
    plot(Median_x,Median_y,'LineWidth',3,'Color',[1.00,0.41,0.16]);
    legend('Average');
end

%% Ideal Array Factor vs BS position Error
z = 1;

R = distance(sensorsPos,BS);

AFIdeal = N; %AF for correct positions

rounds = 500;
%variance1(1,:) = 0.025:0.05:3;
%variance2(1,:) = 3.025:0.05:6;
%variance3(1,:) = 0.025:0.05:0.75;
%variance4(1,:) = 2:0.2:3;
variance = 0.025:0.05:3;
Median_y =  zeros(1,length(variance));
Median_x = zeros(1,length(variance));

for m = 1:1
    AFtest = zeros(rounds,length(variance));
    for i = 1:length(variance)
        for j = 1:rounds*m
            [~, Rerror] = receptor_BSposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,variance(i)*lambda,0);
            AFtest(j,i) = sum(exp(1j*2*pi/lambda*(R-Rerror)));
        end
    end
    AFnormalized = abs(AFtest/AFIdeal);

    figure;
    hd = boxplot(AFnormalized,variance);
    hold on;
    xlabel('Variance (wavelength)');
    ylabel('Normalized Array Factor');
    for i=6:7:7*length(variance)-1
        Median_x(z) = sum(get(hd(i),'XData'))/2;
        Median_y(z) = sum(get(hd(i),'YData'))/2;
        z=z+1;
    end
    plot(Median_x,Median_y,'LineWidth',3,'Color',[1.00,0.41,0.16]);
    legend('Average');
end

%% Drone or BS optimum position 
% Until here drone was located in the center of the field. Now we will be
% located at the optimum position in every simulation

%close all;

[val,BS,R] = optmdrone(fieldx,fieldy,sensorsPos,f,c,N,lambda,25,25);

%x = 0:fieldx/25:fieldx;
%y = fieldy:-fieldy/25:0;
%%
figure;
imagesc(val/max(max(val)))
colorMap = jet(256);
colormap(colorMap);
colorbar;

%% Ideal Signal Amplitude vs BS position Error
z=1;

val = receptor(R,f,c,N,lambda,0);

rounds = 1000;

%variance1(1,:) = 0.025:0.05:3;
%variance2(1,:) = 3.025:0.05:6;
%variance3(1,:) = 0.025:0.05:0.75;
%variance4(1,:) = 2:0.2:3;
variance = 3.025:0.05:4;
Median_y =  zeros(1,length(variance));
Median_x = zeros(1,length(variance));
valtests = zeros(rounds,length(variance));
  
for m = 1:1
    for i = 1:length(variance)
        for j = 1:rounds
            [valtests(j,i),~] = receptor_BSposError(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,variance(i)*lambda,0);
        end
    end


    valnormalized = valtests/val;

    figure;
    hd = boxplot(valnormalized,variance);
    hold on;
    xlabel('Variance (wavelength)');
    ylabel('Normalized received signal amplitude');
    for i=6:7:7*length(variance)-1
        Median_x(z) = sum(get(hd(i),'XData'))/2;
        Median_y(z) = sum(get(hd(i),'YData'))/2;
        z=z+1;
    end
    plot(Median_x,Median_y,'LineWidth',3,'Color',[1.00,0.41,0.16]);
    legend('Average');
end

%% Ideal Array Factor vs BS and sensors' position Error
z=1;

mode = 0;

rounds = 600;

[val,BS,R] = optmdrone(fieldx,fieldy,sensorsPos,f,c,N,lambda,25,25);

%x = 0:fieldx/25:fieldx;
%y = fieldy:-fieldy/25:0;

%figure;
%imagesc(val/max(max(val)))
%colorMap = jet(256);
%colormap(colorMap);
%colorbar;

z = 1;

clear val;

val = receptor(R,f,c,N,lambda,0);


variance = 0.025:0.025:0.5;
Median_y =  zeros(1,length(variance));
Median_x = zeros(1,length(variance));
valtests = zeros(rounds,length(variance));
  
for m = 1:1
    for i = 1:length(variance)
        for j = 1:rounds
            [valtests(j,i),~] = Systemerror(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,variance(i)*lambda,mode);
        end
    end


    valnormalized = valtests/val;

    figure;
    hd = boxplot(valnormalized,variance);
    hold on;
    xlabel('Variance (wavelength)');
    ylabel('Normalized received signal amplitude');
    for i=6:7:7*length(variance)-1
        Median_x(z) = sum(get(hd(i),'XData'))/2;
        Median_y(z) = sum(get(hd(i),'YData'))/2;
        z=z+1;
    end
    plot(Median_x,Median_y,'LineWidth',3,'Color',[1.00,0.41,0.16]);
    legend('Average');
end

%% Ideal Array Factor vs BS and sensors' position Error
mode = 0;

z = 1;

AFIdeal = N;

rounds = 600;

variance = 0.025:0.025:0.5;
Median_y =  zeros(1,length(variance));
Median_x = zeros(1,length(variance));
valtests = zeros(rounds,length(variance));

for m = 1:1
    for i = 1:length(variance)
        for j = 1:rounds
            [~,Rerror] = Systemerror(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,variance(i)*lambda,mode);
            AFtest(j,i) = sum(exp(1j*2*pi/lambda*(R-Rerror)));
        end
    end
    
    AFnormalized = abs(AFtest/AFIdeal);

    figure;
    hd = boxplot(AFnormalized,variance);
    hold on;
    xlabel('Variance (wavelength)');
    ylabel('Normalized Array Factor');
    for i=6:7:7*length(variance)-1
        Median_x(z) = sum(get(hd(i),'XData'))/2;
        Median_y(z) = sum(get(hd(i),'YData'))/2;
        z=z+1;
    end
    plot(Median_x,Median_y,'LineWidth',3,'Color',[1.00,0.41,0.16]);
    legend('Average');
end

%% Signal Adjust
receptorAdjust(R,f,c,N,lambda,1);

%% Ideal Signal Amplitude vs Sensors' position Erro (Amplitude Adjust)
z=1;

R = distance(sensorsPos,BS);

val = receptor(R,f,c,N,lambda,0);
%val = receptorAdjust(R,f,c,N,lambda,0);

rounds = 500;
variance = 0.05:0.025:0.5;
Median_y =  zeros(1,length(variance));
Median_x = zeros(1,length(variance));
valtests = zeros(rounds,length(variance));

for m = 1:1
    for i = 1:length(variance)
        for j = 1:rounds*m
            valtests(j,i) = receptor_SensorposErrorAmpAdjust(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,variance(i)*lambda,0);
        end
    end


    valnormalized = valtests/val;

    figure;
    hd = boxplot(valnormalized,variance);
    hold on;
    xlabel('Variance (wavelength)');
    ylabel('Normalized received signal amplitude');
    for i=6:7:7*length(variance)-1
        Median_x(z) = sum(get(hd(i),'XData'))/2;
        Median_y(z) = sum(get(hd(i),'YData'))/2;
        z=z+1;
    end
    plot(Median_x,Median_y,'LineWidth',3,'Color',[1.00,0.41,0.16]);
    legend('Average');
end

%% Ideal Signal Amplitude vs BS' position Erro (Amplitude Adjust)
z=1;

R = distance(sensorsPos,BS);

%val = receptor(R,f,c,N,lambda,0);
val = receptorAdjust(R,f,c,N,lambda,0);

rounds = 500;
variance = 0.025:0.05:3;
Median_y =  zeros(1,length(variance));
Median_x = zeros(1,length(variance));
valtests = zeros(rounds,length(variance));

for m = 1:1
    for i = 1:length(variance)
        for j = 1:rounds*m
            valtests(j,i) = receptor_BSposErrorAmpAdjust(fieldx,fieldy,BS,R,sensorsPos,f,c,N,lambda,variance(i)*lambda,0);
        end
    end


    valnormalized = valtests/val;

    figure;
    hd = boxplot(valnormalized,variance);
    hold on;
    xlabel('Variance (wavelength)');
    ylabel('Normalized received signal amplitude');
    for i=6:7:7*length(variance)-1
        Median_x(z) = sum(get(hd(i),'XData'))/2;
        Median_y(z) = sum(get(hd(i),'YData'))/2;
        z=z+1;
    end
    plot(Median_x,Median_y,'LineWidth',3,'Color',[1.00,0.41,0.16]);
    legend('Average');
end
