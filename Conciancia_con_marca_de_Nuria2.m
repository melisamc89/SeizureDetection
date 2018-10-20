% Construct binary signals to compare automatic detection with physicians
% detection
clear all
close all
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/26-07-2018/';
directory1='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/26-07-2018/Conciencia/';
directory2='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/Datos/Datos_MATLAB/';
directory3='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/26-07-2018/50s_SlidingWindow/Eigenanalysis/';
cd(directory)
Frequencies=0:0.1:100;
bands=[1,3.5,7.4,12.4,24,97];
Time=1;
patients=['WL';'RF';'MA';'CG';'CM'];
%crisis={[2,3,4,5,6,7,9,11,12],[2,3,4,5,6,7],[1,2,3,4],[37,38,48,50,58,62,66,68],[1,3,6,7,8]};
crisis={[2,3,4,5,6,7],[2,3,4,5,6,7],[1,2,3,4],[37,38,48,50,58,62,66,68],[1,3,6,7,8]};
AOC_value={[2,3,5,4,6,1,0,1,1],[9,9,7,-1,0,-1],[9,9,7,9],[5,5,8,5,5,5,5,-1],[1,-1,6,6,6]};

electrode_onset={{[1,2,3,13,14,15],[4,5,6],[13,14],[14],[3],[2,3,13,14],[3,4],[11,12,13],[11,12,13]},...
    {[25,26,26],[25,26,27],[25,26],[25,26],[26],[25,26]},...
    {[10,11,19,20],[28,29,37,38],[28,29,37,38,46,47],[10,11,19,20]},...
    {[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6]},...
    {[19,10,11],[37,46,47],[37,38,46,47],[19,20,10,11],[19,10,11]}};

electrode_onset_right={{[1,2,3,13,14,15],[4,5,6],[13,14],[14],[3],[2,3,13,14],[3,4],[11,12,13],[11,12,13]},...
    {[],[],[],[],[],[]},...
    {[10,11,19,20],[],[],[10,11,19,20]},...
    {[],[],[],[],[],[],[],[],[]},...
    {[],[37,46,47],[37,38,46,47],[],[]}}; %%12 seizures

electrode_onset_left={{[],[],[],[],[],[],[],[],[]},...
    {[25,26,26],[25,26,27],[25,26],[25,26],[26],[25,26]},...
    {[],[28,29,37,38],[28,29,37,38,46,47],[]},...
    {[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6]},...
    {[19,10,11],[],[],[19,20,10,11],[19,10,11]}}; %% 20 seizures


sampling_rate=[2048,2000,2000,2000,2000];
crisis_time=50;
sampling_rate=[2048,2000,2000,2000,2000];

counter=1;

for p=1:5
    patient_name=patients(p,:);
    name_eigen=strcat(directory3,patient_name,'_all_eigenanalysis');
    load(name_eigen,'-mat')
    last_end=1;
    s_rate=sampling_rate(p);
    for number=1:length(crisis{p})
        day=crisis{p}(number);
        name=strcat(directory2,patient_name,'_',int2str(day));
        load(name,'-mat')
        electrodo_init(number)=SzO/s_rate;
        electrodo_end(number)=SzE/s_rate;
        all_electrodo_init(number)=last_end+SzO/s_rate;
        all_electrodo_end(number)=last_end+SzE/s_rate;
        init_register(number)=last_end;
        last_end=last_end+length(data(1,:))/s_rate;
    end
    init_register(length(crisis{p})+1)=length(all_eigenvalues);
    electrodo_init=floor(electrodo_init)+1;
    electrodo_end=floor(electrodo_end);

    onset_electrodes=electrode_onset{p};
    %onset_electrodes=electrode_onset_right{p};
    %onset_electrodes=electrode_onset_left{p};
    for number=1:length(crisis{p})
        if isempty(onset_electrodes(number))~=1
             start=electrodo_init(number)-100; final=electrodo_init(number)+200;  
            if start>0
                value2=abs(V_crisis_prueba(onset_electrodes{number},:,start:final));
                SzO_mean(counter,:,:)=mean(value2,1);
                AOC(counter)=AOC_value{p}(number);
                counter=counter+1;
            end
        end
    end
end

index=find(AOC>=0);
label=['Delta';'Theta';'Alpha';'Beta1';'Gamma'];
time=-50:1:250;
all_reg_sig=zeros(length(time),5);

for k=1:size(SzO_mean,3)
    %figure(k)
    for i=1:5
        %subplot(2,3,i)
        %scatter(AOC,SzO_mean(:,i,k)','filled')    
        [r m p]=regression(AOC(index),SzO_mean(index,i,k)');
        %hold on
        %x=0:1:10;
        %plot(x,m*x+p,'k','Linewidth',2)
        all_reg(k,i)=r;
        
        for j=1:1000
            [r1 m p]=regression(AOC(index(randperm(length(index)))),SzO_mean(index,i,k)');
            test(j)=r1;
        end
        test=sort(abs(test));
        if abs(r)>test(900)
            %legend(label(i),strcat('P_c_o_e_f_f=',num2str(r),'*'))
            all_reg_sig(k,i)=1;
        %else
            %legend(label(i),strcat('P_c_o_e_f_f=',num2str(r)))
        end
        %xlabel('AOC')
        %ylabel('Mean')
        %xlim([-1 10])
        %ylim([0 1])
        %box on
    end
    %name=strcat(directory1,'Time=',int2str(time(k)));
    %saveas(k,name,'png')
    %close(k)
end

plot(time, all_reg,'Linewidth',2)
hold on
for i=1:size(all_reg,2)
    index=find(all_reg_sig(:,i)>0);
    scatter(time(index),all_reg(index,i));
end


index=find(AOC>=0);
label=['Delta';'Theta';'Alpha';'Beta1';'Gamma'];
time=-50:1:200;
all_reg_sig=zeros(251,5);

[value index2]=sort(AOC(index));
color=colormap(jet(length(index2)));
for i=1:length(index2)
    plot(time,reshape(SzO_mean(index2(i),2,:),[1 251]),'color',color(i,:))
    hold on
end

size_each_patient=[8,4,4,7,4];
last=1;
used_AOC=AOC(index);
color=colormap(jet(11));
for p=1:5
    subplot(2,3,p)
    group=used_AOC(last:last+size_each_patient(p)-1);
    [new_group index2]=sort(group);
    for i=1:size_each_patient(p)
        plot(time,reshape(SzO_mean(last+index2(i)-1,5,:),[1 251]),'color',color(new_group(index2(i))+1,:))
        hold on
    end
    last=last+size_each_patient(p);
end
