function [t, adjustsignals, adjustsignals_Loss, receivedsignal,receivedsignal_Loss, PropagationLoss] = receptorgroups(group,R,Rmax,f,c,N,lambda,mode)
    
    ngroups = max(group(:,4));
    Num = zeros(1,ngroups);
    for i = 1:ngroups
        Num(i) = size(group(group(:,4)==i),1);
    end

    %To remove z component error from sensors -> activate line below
    phase = (Rmax-R).*2*pi/lambda;
    traveling_time = R/c;
    %size(traveling_time)
    %The receiver will wait a time interval equal to t for all signals
    tmax = max(max(traveling_time))*3;
    t=linspace(0,tmax,1000000); % 1 período
    %tstem = 0:1/f:tmax;

    %signals = ones(N,length(t));
    signals = zeros(100,length(t));
    % for i = 1:N
    %     if mode == 1
    %         signals(i,:) = cos(2*pi*f*t-phase(i));
    %     elseif mode == 0
    %         signals(i,:) = cos(2*pi*f*t);
    %     end
    %     %AF = AF + exp(1j*2*pi*(-R(i))/lambda);
    % end
    
    for i = 1:ngroups
        for j = 1:Num(i)
            if mode == 1
                signals(j+20*(i-1),:) = cos(2*pi*f*t-phase(j+20*(i-1)));
            elseif mode == 0
                signals(j,:) = cos(2*pi*f*t);
            end
            %AF = AF + exp(1j*2*pi*(-R(i))/lambda);
        end
    end
    
    PropagationLoss = (4*pi*R/lambda).^2;
    
    [adjustsignals, adjustsignals_Loss, receivedsignal,receivedsignal_Loss] = adjsig(size(group,1),t,signals,traveling_time,PropagationLoss);
    
    % Signals plot without time adjust
    %figure(10);
    %plot(t,signals(:,:));