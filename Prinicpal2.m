p=1;
day=3;
load('wagner_3','-mat')
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/EpileptogenicityIndex/'
s_rate=header.samplingrate;
Fs = s_rate;            % Sampling frequency
c = colormap(jet(44));
channels=ChannelNumber(patients{p},day,header)
for i=1:44
    signal=data(i,:);
    window_size = s_rate*1;
    overlap = floor(s_rate*0.5);
    [Y,F,T,PSD] = spectrogram(signal,window_size,overlap,0:0.1:100,s_rate,'power','yaxis');
    TT = linspace(0,length(signal)/s_rate,length(signal));
    
    figure(1)
    EI=EpileptogenicityIndex(F,0.1,T,0.5,PSD);
    plot(T,EI+i)
    hold on
   
    figure(2)
    CEI=CumulativeEpileptogenicityIndex(EI,0);
    plot(T,CEI+i)
    hold on
    
end
figure(1)
 xlim([0 time(end)])
 line([patients{p}.crisis.electrodo_init(day)/s_rate patients{p}.crisis.electrodo_init(day)/s_rate],...
        [0 80],'Color','k','Linewidth',1)
 line([patients{p}.crisis.propagation(day)/s_rate patients{p}.crisis.propagation(day)/s_rate],...
         [0 80],'Color','g','Linewidth',1)
 line([patients{p}.crisis.AOC(day)/s_rate patients{p}.crisis.AOC(day)/s_rate],...
         [0 80],'Color','r','Linewidth',1)
 line([patients{p}.crisis.electrodo_end(day)/s_rate patients{p}.crisis.electrodo_end(day)/s_rate],...
        [0 80],'Color','k','Linewidth',1)
 xlabel('Time [s]');
 ylabel('EIndex');
    
 figure(2)
 ylimits=[min(CEI) max(CEI)+i];
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
    
    saveas(1,strcat(directory,'EI_Channels_crisis6','png')
    close(1)
    saveas(2,strcat(directory,'EI_Channels_crisis5','png')
    close(2)