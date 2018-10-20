% Construct binary signals to compare automatic detection with physicians
% detection
clear all
close all
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/26-07-2018/';
directory2='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/Datos/Datos_MATLAB/';
cd(directory)
Frequencies=0:0.1:100;
bands=[1,3.5,7.4,12.4,24,97];
Time=1;
patients=['WL';'RF';'MA';'CG';'CM'];
crisis={[2,3,4,5,6,7,9,11,12],[2,3,4,5,6,7],[1,2,3,4],[37,38,48,50,58,62,66,68,68],[1,3,6,7,8]};
electrode_onset={{[1,2,3,13,14,15],[4,5,6],[13,14],[14],[3],[2,3,13,14],[3,4],[11,12,13],[11,12,13]},...
    {[25,26,26],[25,26,27],[25,26],[25,26],[26],[25,26]},...
    {[10,11,19,20],[28,29,37,38],[28,29,37,38,46,47],[10,11,19,20]},...
    {[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6]},...
    {[1,2,3,19,10,11],[37,46,47],[28,29,37,38,46,47],[19,20,10,11],[19,10,11]}};
sampling_rate=[2048,2000,2000,2000,2000];
crisis_time=25;
sampling_rate=[2048,2000,2000,2000,2000];

overlapping=0;
counter=0;
false_possitive=0;
false_negative=0;
positive_negative=0;
PRES=900;
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
    
    for number=1:length(crisis{p})
        onset_electrodes=electrode_onset{p};
        binary_signal_2{crisis{p}(number)}=zeros(length(onset_electrodes{number}),length(init_register(number):init_register(number+1)));
        binary_signal{crisis{p}(number)}=zeros(1,length(init_register(number):init_register(number+1)));
        binary_signal{crisis{p}(number)}(1,electrodo_init(number):electrodo_end(number))=1;
        for j=1:length(onset_electrodes{number})
             sort_value=sort(test_distribution(onset_electrodes{number}(j),:));
             test_value=sort_value(PRES);
             index=find(all_eigenvalues(onset_electrodes{number}(j),floor(init_register(number)):floor(init_register(number+1)))>test_value);
             binary_signal_2{crisis{p}(number)}(j,index)=1;
        end
    end
    clear i
    for k=1:length(binary_signal)
        for j=1:size(binary_signal_2{k},1)
            if isempty(binary_signal_2{k})~=1 && isempty(binary_signal{k})~=1 
            overlapping=overlapping+length(find(binary_signal{k}-binary_signal_2{k}(j,:)));
            counter=counter+length(binary_signal{k});
            differences=find(binary_signal{k}~=binary_signal_2{k}(j,:));
            false_possitive=false_possitive+length(find(binary_signal_2{k}(j,differences)==1));
            false_negative=false_negative+length(find(binary_signal_2{k}(j,differences)==0));
            positive_negative=positive_negative+length(find(binary_signal{k}(differences)==0));
            end
        end
    end
    clear binary_signal binary_signal_2 electrodo_init electrodo_end init_register last_end
end
performance=abs(overlapping)/counter;
false_possitive=false_possitive/counter;
false_negative=false_negative/counter;
        