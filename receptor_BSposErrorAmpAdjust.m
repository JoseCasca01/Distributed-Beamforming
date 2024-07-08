function [val,Rerror] = receptor_BSposErrorAmpAdjust(fieldx,fieldy,BS,R,group,f,c,N,lambda,erro,mode)

    %%The receiver will wait a time interval equal to t for all signals
    t=linspace(0,max(R)*3/c,10000); % 1 period
    
    PT = 1;
    GT = 1;
    GR = 1;
    
    PR = zeros(N,1);
    for i = 1:N
        %Friis formula for free-space propagation    
        PR(i) = PT * GT * GR * lambda^2/(16*pi^2*R(i)^2);
    end
    
    %The sensor furthest from the receiver will be the reference sensor with phase 0
    phaseCorrect = (max(R)-R).*2*pi/lambda;
    traveling_time = R/c;

    PropagationLoss = (4*pi*R/lambda).^2;
    BSerror = BS + (erro*randn(1,1));
    
    Rerror = distance(group,BSerror);    
    SignalsNewAmp = AmpAdjust(R);
    SignalsNewAmperror = AmpAdjust(Rerror);

    phase = (max(Rerror)-Rerror).*2*pi/lambda;
    
    signals = zeros(N,length(t));
    signalsCorrect = zeros(N,length(t));

    for i = 1:N
        signals(i,:) = SignalsNewAmperror(i)*cos(2*pi*f*t-phase(i));
        signalsCorrect(i,:) = SignalsNewAmp(i)*cos(2*pi*f*t-phaseCorrect(i));
    end

    [~, adjustsignals_Loss, receivedsignal, receivedsignal_Loss] = adjsig(N,t,signals,traveling_time,PropagationLoss);
    [~, ~, ~, idealsignal_Loss] = adjsig(N,t,signalsCorrect,traveling_time,PropagationLoss);

    aux = sum(adjustsignals_Loss);
    val = max(aux(t>=1.1*traveling_time(R==max(R))));
    %disp(aux);
    %disp(val);
    if mode == 1
        close all;

        figure(1);
        %%%%% 2D%%%%%
        plot(BS(1),BS(2),'X'),hold on;
        plot(group(:,1),group(:,2),'O')
        plot(BSerror(:,1),BSerror(:,2),'diamond');
        title('Field');
        ylabel('yfield (m)');
        xlabel('xfield (m)');
        legend('Drone', 'Node Position', 'Drone Position with error');
        axis([0, fieldy, 0, fieldx]);
        %%%%%3D%%%%%
        figure(8);
        plot3(BS(1),BS(2),BS(3),'X'),hold on;
        plot3(group(:,1),group(:,2),group(:,3),'O')
        plot3(BSerror(:,1),BSerror(:,2),BSerror(:,3),'diamond');
        title('Field');
        ylabel('yfield (m)');
        xlabel('xfield (m)');
        legend('Drone', 'Node Position', 'Drone Position with error');
        axis([0, fieldy, 0, fieldx]);

        figure(2);
        plot(1:N,PropagationLoss,'-X');
        title('Propagation Loss per sensor')
        ylabel('P_R');
        xlabel('Sensor number');
    
        % figure(3);
        % plot(t,adjustsignals(:,:));
        % title('Signals Received With Position Error');
        % xlabel('t(s)');
        % ylabel('Amplitude');
        
        figure(4);
        plot(t,receivedsignal);
        title('Received Signal');
        xlabel('t(s)');
        ylabel('Amplitude');
        axis([0 max(t)*1.01 min(receivedsignal)-1 max(receivedsignal)+1]);
        
        % figure(5);
        % plot(t,adjustsignals_Loss);
        % title('Signals Received with Losses With Position Error');
        % xlabel('t(s)');
        % ylabel('Amplitude');
        
        figure(6);
        plot(t,receivedsignal_Loss);
        title('Received Signal with Losses');
        xlabel('t(s)');
        ylabel('Amplitude');
        axis([0 max(t)*1.01 min(receivedsignal_Loss)*1.01 max(receivedsignal_Loss)*1.01]);

        figure(7);
        plot(t,idealsignal_Loss-receivedsignal_Loss);
        title('Received Signals Sum Differences');
        xlabel('t(s)');
        ylabel('Amplitude');
        
        %Signals plot without time adjust
        %figure(10);
        %plot(t,signals(:,:));

        disp(val)
    end
end