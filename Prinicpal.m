p=1;
day=3;
name=strcat('wagner_',int2str(day));
load(name,'-mat')
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/EpileptogenicityIndex/'
s_rate=header.samplingrate;
Fs = s_rate;            % Sampling frequency
c = colormap(jet(44))
channels=ChannelNumber(patients{p},day,header)
for i=1:44
    signal=data(i,:);
    figure(i)
    subplot(4,1,1)
    if find(channels==i)
        title('INIT CHANNEL')
    end
    set(gca,'fontsize',8)
    hold on
    time=(0:length(data(i,:))-1)/s_rate;
    plot(time,signal,'b','Linewidth',1)
    xlim([0 time(end)])
    ylim([min(signal) max(signal)])
    hold on
    line([patients{p}.crisis.electrodo_init(day)/s_rate patients{p}.crisis.electrodo_init(day)/s_rate],...
        [min(signal) max(signal)],'Color','k','Linewidth',3)
    line([patients{p}.crisis.propagation(day)/s_rate patients{p}.crisis.propagation(day)/s_rate],...
        [min(signal) max(signal)],'Color','g','Linewidth',4)
    line([patients{p}.crisis.AOC(day)/s_rate patients{p}.crisis.AOC(day)/s_rate],...
        [min(signal) max(signal)],'Color','r','Linewidth',2)
    line([patients{p}.crisis.electrodo_end(day)/s_rate patients{p}.crisis.electrodo_end(day)/s_rate],...
        [min(signal) max(signal)],'Color','k','Linewidth',3)
    ylabel('LFP [uV]')
    
    subplot(4,1,2)
    set(gca,'fontsize',8)
    hold on
    window_size = s_rate*1;
    overlap = floor(s_rate*0.5);
    [Y,F,T,PSD] = spectrogram(signal,window_size,overlap,0:0.1:100,s_rate,'yaxis');
    TT = linspace(0,length(signal)/s_rate,length(signal));
    imagesc(TT,F,(abs(PSD).^0.25))
    set(gca,'YDir','normal'), colormap(jet)
    axis xy
    axis tight
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    %set(gca,'XTick',[0 1200 2400 3600],...
    %    'YTick',[0 2 4 6 8 10 12 14 16 18])
    box on
    hold on
    line([patients{p}.crisis.electrodo_init(day)/s_rate patients{p}.crisis.electrodo_init(day)/s_rate],...
        [min(F) max(F)],'Color','k','Linewidth',3)
    line([patients{p}.crisis.propagation(day)/s_rate patients{p}.crisis.propagation(day)/s_rate],...
        [min(F) max(F)],'Color','g','Linewidth',4)
    line([patients{p}.crisis.AOC(day)/s_rate patients{p}.crisis.AOC(day)/s_rate],...
        [min(F) max(F)],'Color','r','Linewidth',2)
    line([patients{p}.crisis.electrodo_end(day)/s_rate patients{p}.crisis.electrodo_end(day)/s_rate],...
        [min(F) max(F)],'Color','k','Linewidth',3)
    legend('Init','Propagation','AOC','End')
    
    subplot(4,1,3)
    EI=EpileptogenicityIndex(F,0.1,T,0.5,PSD);
    plot(T,EI)
    hold on
    ylim([0 4])
    xlim([0 time(end)])
    line([patients{p}.crisis.electrodo_init(day)/s_rate patients{p}.crisis.electrodo_init(day)/s_rate],...
        [0 6],'Color','k','Linewidth',1)
    line([patients{p}.crisis.propagation(day)/s_rate patients{p}.crisis.propagation(day)/s_rate],...
         [0 6],'Color','g','Linewidth',1)
    line([patients{p}.crisis.AOC(day)/s_rate patients{p}.crisis.AOC(day)/s_rate],...
         [0 6],'Color','r','Linewidth',1)
    line([patients{p}.crisis.electrodo_end(day)/s_rate patients{p}.crisis.electrodo_end(day)/s_rate],...
        [0 6],'Color','k','Linewidth',1)
    xlabel('Time [s]');
    ylabel('EIndex');
  
    subplot(4,1,4)
    CEI=CumulativeEpileptogenicityIndex(EI,0);
    plot(T,CEI)
    hold on
    ylimits=[min(CEI) max(CEI)];
    ylim(ylimits)
    xlim([0 time(end)])
    line([patients{p}.crisis.electrodo_init(day)/s_rate patients{p}.crisis.electrodo_init(day)/s_rate],...
        ylimits,'Color','k','Linewidth',1)
    line([patients{p}.crisis.propagation(day)/s_rate patients{p}.crisis.propagation(day)/s_rate],...
         ylimits,'Color','g','Linewidth',1)
    line([patients{p}.crisis.AOC(day)/s_rate patients{p}.crisis.AOC(day)/s_rate],...
         ylimits,'Color','r','Linewidth',1)
    line([patients{p}.crisis.electrodo_end(day)/s_rate patients{p}.crisis.electrodo_end(day)/s_rate],...
        ylimits,'Color','k','Linewidth',1)
    xlabel('Time [s]');
    ylabel('CumEI');
    
    
    saveas(i,strcat(directory,'EI_Channel_',header.channels(i,5:end)),'png')
    close(i)
end