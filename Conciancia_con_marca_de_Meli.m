% Construct binary signals to compare automatic detection with physicians
% detection
clear all
close all
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/26-07-2018/';
directory1='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/26-07-2018/50s_SlidingWindow/Eigenanalysis/';
directory2='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/Datos/Datos_MATLAB/';
cd(directory)
Frequencies=0:0.1:100;
bands=[1,3.5,7.4,12.4,24,97];
Time=1;
patients=['WL';'RF';'MA';'CG';'CM'];
crisis={[2,3,4,5,6,7,9,11,12],[2,3,4,5,6,7],[1,2,3,4],[37,38,48,50,58,62,66,68],[1,3,6,7,8]};
AOC_value={[0,3,5,4,6,1,0,1,1],[9,9,7,-1,0,-1],[9,9,7,9],[5,5,8,5,5,5,5,-1],[1,-1,6,6,6]};
AOC=[0,3,5,4,6,1,0,1,1,9,9,7,-1,0,-1,9,9,7,9,5,5,8,5,5,5,5,-1,1,-1,6,6,6];
electrode_onset={{[1,2,3,13,14,15],[4,5,6],[13,14],[14],[3],[2,3,13,14],[3,4],[11,12,13],[11,12,13]},...
    {[25,26,26],[25,26,27],[25,26],[25,26],[26],[25,26]},...
    {[10,11,19,20],[28,29,37,38],[28,29,37,38,46,47],[10,11,19,20]},...
    {[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6]},...
    {[19,10,11],[37,46,47],[37,38,46,47],[19,20,10,11],[19,10,11]}};

electrode_onset={{[1:44],[1:44],[1:44],[1:44],[1:44],[1:44],[1:44],[1:44],[1:44]},...
   {[1:45],[1:45],[1:45],[1:45],[1:45],[1:45]},...
    {[1:54],[1:54],[1:54],[1:54]},...
    {[1:39],[1:39],[1:39],[1:39],[1:39],[1:39],[1:39],[1:39],[1:39]},...
    {[1:54],[1:54],[1:54],[1:54],[1:54]}};


sampling_rate=[2048,2000,2000,2000,2000];
crisis_time=25;
sampling_rate=[2048,2000,2000,2000,2000];

counter=1;

for p=1:5
    patient_name=patients(p,:);
    name_eigen=strcat(directory,patient_name,'_all_eigenanalysis');
    load(name_eigen,'-mat')
    last_end=1;
    s_rate=sampling_rate(p);
    for number=1:length(crisis{p})
        day=crisis{p}(number);
        name=strcat(directory2,patient_name,'_',int2str(day));
        load(name,'-mat')
        if p==1 && number==1 
            save_channels=channels;end
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
    if p==1 
        channels=save_channels(:,5:8);end
    
    N=length(channels);
    if p==1
        N=N-2; end
    for i=1:N
        sorted_test=sort(test_distribution(i,:));
        baseline=sorted_test(900);
        for number=1:length(crisis{p})
            recording=all_eigenvalues(i,floor(init_register(number))+1:floor(init_register(number+1)));
            mark1=find(recording>baseline);
            in_line=[electrodo_init(number)-50:1:electrodo_init(number)+200];
            [tf,loc]=ismember(in_line,mark1);
            if sum(tf)~=0
                mark=min(in_line(find(tf>0,1)))+floor(init_register(number))-50;
            else
                mark=0;
            end
            index(i,number)=mark;
        end
    end
    onset_electrodes=electrode_onset{p};
    for number=1:length(crisis{p})
        value1=0; value2=0;value3=0;
        for i=1:length(onset_electrodes{number})
            if index(i,number)
            start=index(i,number); final=index(i,number)+200;                 
            before_start=start-10; 
            if before_start<0
               before_start=1;end
            if final>length(all_eigenvalues)
                final=length(all_eigenvalues);end
            value1=value1+mean(abs(V_crisis_prueba(onset_electrodes{number}(i),:,start:final)),3);
            value2=value2+mean(abs(V_crisis_prueba(onset_electrodes{number}(i),:,before_start:start)),3);
            value3=value3+mean(abs(V_crisis_prueba(onset_electrodes{number}(i),:,start:start+10)),3);
            end         
        end
        all_Seizure_mean(counter,:)=value1/length(onset_electrodes{number});
        before_SzO_mean(counter,:)=value2/length(onset_electrodes{number});
        after_SzO_mean(counter,:)=value3/length(onset_electrodes{number});
        AOC(counter)=AOC_value{p}(number);
        counter=counter+1;
    end
end

index=find(AOC>0);
label=['Delta';'Theta';'Alpha';'Beta1';'Gamma'];
for i=1:5
    subplot(3,5,i)
    scatter(AOC,all_Seizure_mean(:,i),'filled')    
    [r m p]=regression(AOC(index),all_Seizure_mean(index,i)');
    hold on
    x=0:1:10;
    plot(x,m*x+p,'k','Linewidth',2)
    
    for j=1:1000
        [r1 m p]=regression(AOC(index(randperm(length(index)))),all_Seizure_mean(index,i)');
        test(j)=r1;
    end
    test=sort(abs(test));
    if abs(r)>test(900)
        legend(label(i),strcat('P_c_o_e_f_f=',num2str(r),'*'))
    else
        legend(label(i),strcat('P_c_o_e_f_f=',num2str(r)))
    end
    xlabel('AOC')
    ylabel('Mean')
    xlim([-1 10])
    ylim([0 1])
    box on
    
    subplot(3,5,i+5)
    scatter(AOC,before_SzO_mean(:,i),'filled')
    [r m p]=regression(AOC(index),before_SzO_mean(index,i)');
    hold on
    x=0:1:10;
    plot(x,m*x+p,'k','Linewidth',2)
    for j=1:1000
        [r1 m p]=regression(AOC(index(randperm(length(index)))),before_SzO_mean(index,i)');
        test(j)=r1;
    end
    test=sort(abs(test));
    if abs(r)>test(900)
        legend(label(i),strcat('P_c_o_e_f_f=',num2str(r),'*'))
    else
        legend(label(i),strcat('P_c_o_e_f_f=',num2str(r)))
    end
    xlabel('AOC')
    ylabel('Mean')
    xlim([-1 10])
    ylim([0 1])
    box on
    
    subplot(3,5,i+10)
    scatter(AOC,after_SzO_mean(:,i),'filled')
    [r m p]=regression(AOC(index),after_SzO_mean(index,i)');
    hold on
    x=0:1:10;
    plot(x,m*x+p,'k','Linewidth',2)
    for j=1:1000
        [r1 m p]=regression(AOC(index(randperm(length(index)))),after_SzO_mean(index,i)');
        test(j)=r1;
    end
    test=sort(abs(test));
    if abs(r)>test(900)
        legend(label(i),strcat('P_c_o_e_f_f=',num2str(r),'*'))
    else
        legend(label(i),strcat('P_c_o_e_f_f=',num2str(r)))
    end
    xlabel('AOC')
    ylabel('Mean')
    xlim([-1 10])
    ylim([0 1])
    box on
end
    