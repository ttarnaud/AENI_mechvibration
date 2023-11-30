function [RMS,RMSabs, Q2, Q2abs] = calcRMS(Signal,InputSignal,timeReconS,mirror,disp)
%calculate root mean square between reconstructed and input signal

InS = interp1(InputSignal(1,:),InputSignal(2,:),timeReconS);
% Check if signal lengths are the same
if length(Signal)~=length(InS)
    error('signal lengths not the same')
end
Signal = (-1)^mirror.*Signal;

% redefined normalization in such way that sign diffrence doesn't
% affect RMS
%Signal_peakval = max(abs(Signal));
Signal_peakval = Signal(max(abs(Signal))==abs(Signal));
Signal_peakval = Signal_peakval(1);
%InS_peakval = max(abs(InS));
InS_peakval = InS(max(abs(InS))==abs(InS));
InS_peakval = InS_peakval(1);

% RMS
RMS = sqrt(mean((Signal/Signal_peakval-InS/InS_peakval).^2));
RMS2 = sqrt(mean((Signal/(-Signal_peakval)-InS/InS_peakval).^2));
if RMS2<RMS
    %by chance (due to low sampling) in a periodic signal (eg sin)
    %the negative peak could have a higher value than positive peak
    %but swap not necessary
    RMS = RMS2;
end
%RMSabs
RMSabs = sqrt(mean((abs(Signal/Signal_peakval)-abs(InS/InS_peakval)).^2));

% Q2
Q2 = sqrt(mean((Signal/Signal_peakval-InS/InS_peakval).^2))./sqrt(sum((InS/InS_peakval).^2));
Q22 = sqrt(mean((Signal/(-Signal_peakval)-InS/InS_peakval).^2))./sqrt(sum((InS/InS_peakval).^2));
if Q22<Q2
    %explanation see RMS
    Q2 = Q22;
end
%Q2abs
Q2abs = sqrt(mean((abs(Signal/Signal_peakval)-abs(InS/InS_peakval)).^2))./sqrt(sum((abs(InS/InS_peakval)).^2));

if disp
    figure
    subplot(2,1,1)
    plot(timeReconS,(Signal/max(abs(Signal))),'DisplayName','CompleteSignal')
    hold on
    plot(timeReconS,InS/max(abs(InS)),'DisplayName','SampledInput')
    hold off
    ylabel('normalized signals')
    legend('show')
    subplot(2,1,2)
    plot(timeReconS,(abs(Signal)/max(abs(Signal))),'DisplayName','CompleteSignal')
    hold on
    plot(timeReconS,abs(InS)/max(abs(InS)),'DisplayName','Sampled Input')
    hold off
    ylabel('normalized & rectified signals')
    xlabel('time [ms]')
    legend('show')
    
end
end