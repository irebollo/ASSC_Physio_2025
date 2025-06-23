function [figure_EGGtimeseries] = plot_EGGtimeseries(EGG_raw)

%%% Edited for Brain-Body Waves 2024 LB 

figure_EGGtimeseries = figure
nChannels = length(EGG_raw.label);
for iChannel = 1:nChannels
    subplot(nChannels,1,iChannel)
    plot(EGG_raw.time{1,1}(1,:),EGG_raw.trial{1,1}(iChannel,:))
    title(strcat(('EGG timeseries n'),num2str(iChannel)),'fontsize',11);
    xlabel('time in s')
    ylabel('amplitude')
    set(gca,'fontsize',12)
    set(gcf,'Position',[1000 1000 1500 15000])
end

end