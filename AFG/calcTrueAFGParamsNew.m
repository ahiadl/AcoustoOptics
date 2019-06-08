function [AFG] = generateSignalsForAFG(vars)
% this function generate acousto optic external clk and a single pulse of
% US buily of given number of sine function cycles and padded with zeros.

% Generate The Signal
dt = 1/vars.fSclk;
t = (0:1:(vars.SclkSamplesInTrain-1)) * dt;

pulseData = zeros(samplesInTrain, 1);
pulseData(1:SclkSamplesInPulse) = sin(2*pi*fSin*t(1:SclkSamplesInPulse));

if (AFG.debug)
    figure()
    stem(pulseData)
    title('Signal Downloaded to AFG');
end

%Generate The Clock
clkData      = ones(1,extClkSamplesPerCyc);
dutyCycleIdx = floor(extClkSamplesCycle*(AFG.extClkDcyc/100))+1;
clkData(dutyCycleIdx:end) = 0;
clkData = repmat(clkData, 1, extClkCycles);

if (AFG.debug)
    figure()
    stem(clkData)
    title('Signal Downloaded to AFG');
    
    figure(); 
    plot(linspace(0,AFG.tTrainTrue,samplesInTrain)*1e6, pulseData)
    set(gca, 'FontSize', 20)
    title('Input US pulse')
    xlabel('t[\mus]')
    ylabel('Amplitude')
    xlim([0, 5])
end

AFG.pulseData = pulseData;
AFG.clkData   = clkData;

end
    
   