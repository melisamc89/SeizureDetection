%Construct binary signals to compare automatic detection with physicians
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
%crisis={[2,3,4,5,6,7,9,11,12],[2,3,4,5,6,7],[1,2,3,4],[37,38,48,50,58,62,66,68],[1,3,6,7,8]};
crisis={[2,3,4,5,6,7],[2,3,4,5,6,7],[1,2,3,4],[37,38,48,50,58,62,66,68],[1,3,6,7,8]};

AOC_value={[0,3,5,4,6,1,0,1,1],[9,9,7,-1,0,-1],[9,9,7,9],[5,5,8,5,5,5,5,-1],[1,-1,6,6,6]};

%AOC=[0,3,5,4,6,1,0,1,1,9,9,7,-1,0,-1,9,9,7,9,5,5,8,5,5,5,5,-1,1,-1,6,6,6];
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
        electrodo_init(number)=SzO/s_rate-10;
        electrodo_end(number)=SzE/s_rate+10;
        all_electrodo_init(number)=last_end+SzO/s_rate;
        all_electrodo_end(number)=last_end+SzE/s_rate;
        init_register(number)=last_end;
        last_end=last_end+length(data(1,:))/s_rate;
    end
    init_register(length(crisis{p})+1)=length(all_eigenvalues);
    electrodo_init=max(floor(electrodo_init)+1,1);
    electrodo_end=max(floor(electrodo_end),1);  
    init_register(length(crisis{p})+1)=length(all_eigenvalues);
    electrodo_init=max(floor(electrodo_init)+1,1);
    electrodo_end=max(floor(electrodo_end),1);    
    for number=1:length(crisis{p})
        onset_electrodes=electrode_onset{p};
        binary_signal_2{number}=zeros(length(onset_electrodes{number}),length(init_register(number):init_register(number+1)));
        binary_signal{number}=zeros(1,length(init_register(number):init_register(number+1)));
        binary_signal{number}(1,electrodo_init(number):electrodo_end(number))=1;
        for j=1:length(onset_electrodes{number})
             sort_value=sort(test_distribution(onset_electrodes{number}(j),:));
             test_value=sort_value(950);
             index=find(all_eigenvalues(onset_electrodes{number}(j),floor(init_register(number)):floor(init_register(number+1)))>test_value);
             binary_signal_2{number}(j,index)=1;
        end
    end
    
    for k=1:length(binary_signal)
        for j=1:size(binary_signal_2{k},1)
            if isempty(binary_signal_2{k})~=1 && isempty(binary_signal{k})~=1       
                s_real=find(binary_signal{k});
                s_detected=find(binary_signal_2{k}(j,:));
                [tf,loc]=ismember(s_real,s_detected);
                if length(find(tf))>crisis_time+crisis_time*1 || length(find(tf))==length(s_real)
                    mark=min(s_real(find(tf)))+init_register(k);
                else
                    mark=0;
                end
            end 
            index(j,k)=mark;
        end
    end 
    onset_electrodes=electrode_onset{p};
    %onset_electrodes=electrode_onset_right{p};
    %onset_electrodes=electrode_onset_left{p};
    for number=1:length(crisis{p})
        if isempty(onset_electrodes(number))~=1
            counter2=1;
            detection=0;
            for j=1:length(onset_electrodes{number})
                start=index(j,number)-100; final=index(j,number)+200;
                if start>0
                    value2(counter2,:)=abs(all_eigenvalues(onset_electrodes{number}(j),start:final));
                    counter2=counter2+1;
                    detection=1;
                end
            end
            if detection==1
            SzO_mean(counter,:)=mean(value2,1);
            AOC(counter)=AOC_value{p}(number);
            counter=counter+1;
            end
        end
    end
    clear binary_signal binary_signal_2
end


time=-75:1:225;
figure(2)
[new_AOC new_index]=sort(AOC(find(AOC>=0)));
color=colormap(jet(length(new_AOC)));
for i=1:length(new_AOC)
    plot(time,SzO_mean(new_AOC(i),:),'color',color(i,:))
    hold on
end

time=-75:1:225;
all_reg_sig=zeros(301,5);
index=find(AOC>=0);
for k=1:size(SzO_mean,2)
   [r m p]=regression(AOC(index),SzO_mean(index,k)');
   all_reg(k)=r;
   for j=1:1000
      [r1 m p]=regression(AOC(index(randperm(length(index)))),SzO_mean(index,k)');
      test(j)=r1;
   end
  test=sort(abs(test));
  if abs(r)>test(950)
      all_reg_sig(k)=1;
  end
end

figure(2)
plot(time, all_reg,'Linewidth',2)
hold on
index=find(all_reg_sig(:)>0);
scatter(time(index),all_reg(index));

%%

index=find(AOC>=0);
%label=['Delta';'Theta';'Alpha';'Beta1';'Gamma'];
time=-75:1:225;
all_reg_sig=zeros(length(time),5);

for k=70:150
    figure(k)
    scatter(AOC,SzO_mean(:,k)','filled')    
    [r m p]=regression(AOC(index),SzO_mean(index,k)');
    hold on
    x=0:1:10;
    plot(x,m*x+p,'k','Linewidth',2)
    all_reg(k)=r;
    for j=1:1000
       [r1 m p]=regression(AOC(index(randperm(length(index)))),SzO_mean(index,k)');
       test(j)=r1;
    end
    test=sort(abs(test));
    if abs(r)>test(900)
        legend(strcat('P_c_o_e_f_f=',num2str(r),'*'))
            all_reg_sig(k)=1;
        else
            legend(strcat('P_c_o_e_f_f=',num2str(r)))
    end
    xlabel('AOC')
    ylabel('Eigenvalue')
    xlim([-1 10])
    ylim([0 100])
    box on
    name=strcat(directory,'Time=',int2str(time(k)),'s');
    saveas(k,name,'png')
    close(k)
end