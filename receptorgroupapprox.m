function [t, adjustsignals, adjustsignals_Loss, receivedsignal,receivedsignal_Loss, PropagationLoss] = receptorgroupapprox(group,BS,f,c,N,lambda,mode)
    
    ngroups = max(group(:,4));
    masscenter = zeros(ngroups,3);
    Num = zeros(1,ngroups);

    for i = 1:ngroups
        masscenter(i,1) = sum(group(group(:,4)==i,1)) / size(group(group(:,4)==i),1);
        masscenter(i,2) = sum(group(group(:,4)==i,2)) / size(group(group(:,4)==i),1);
        Num(i) = size(group(group(:,4)==i),1);
    end

    masscenter(:,3) = 0;
    %size(masscenter)
    masscenter
    R = distance(masscenter,BS);
    %size(R)

    phase = (max(R)-R).*2*pi/lambda;
    %size(phase)
    phase
    Rsensors = distance(group(:,1:3),BS);

    traveling_time = Rsensors/c;
    %traveling_time
    %size(traveling_time)
    %max(traveling_time)
    %The receiver will wait a time interval equal to t for all signals
    tmax = max(traveling_time)*3;
    t=linspace(0,tmax,1000000); % 1 período
    %tstem = 0:1/f:tmax;

    %signals = ones(N,length(t));
    signals = zeros(sum(Num),length(t));

    for i = 1:ngroups
        for j = 1:Num(i)
            if mode == 1
                signals(j+20*(i-1),:) = cos(2*pi*f*t-phase(i));
            elseif mode == 0
                signals(j,:) = cos(2*pi*f*t);
            end
            %AF = AF + exp(1j*2*pi*(-R(i))/lambda);
        end
    end
    
    %size(signals)
    PropagationLoss = (4*pi*Rsensors/lambda).^2;
    
    %adjustsignals = zeros(ngroups,length(t));
    %adjustsignals_Loss = zeros(ngroups,length(t)); 
    %receivedsignal = zeros(ngroups,length(t));
    %receivedsignal_Loss = zeros(ngroups,length(t));

    [adjustsignals, adjustsignals_Loss, receivedsignal,receivedsignal_Loss] = adjsig(size(group,1),t,signals,traveling_time,PropagationLoss);
    %size(receivedsignal_Loss)

end