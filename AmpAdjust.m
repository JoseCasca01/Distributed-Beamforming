function [signalsNewAmp] = AmpAdjust(R)
    
    aux = max(R);
    signalsNewAmp = zeros(1,length(R));
    signalsNewAmp(R==aux) = 1;

    for i = 1:length(R)
        signalsNewAmp(i) = (R(i)/aux)^2;
    end
end