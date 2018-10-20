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
    {[19,10,11],[37,46,47],[37,38,46,47],[19,20,10,11],[19,10,11]}};

electrode_onset={{[1:44],[1:44],[1:44],[1:44],[1:44],[1:44],[1:44],[1:44],[1:44]},...
    {[1:45],[1:45],[1:45],[1:45],[1:45],[1:45]},...
    {[1:54],[1:54],[1:54],[1:54]},...
    {[1:39],[1:39],[1:39],[1:39],[1:39],[1:39],[1:39],[1:39],[1:39]},...
    {[1:54],[1:54],[1:54],[1:54],[1:54]}};

sampling_rate=[2048,2000,2000,2000,2000];
crisis_time=50;
sampling_rate=[2048,2000,2000,2000,2000];
total_seizure=33;

detected_seizure=0;
predicted_seizure=0;
counter=0;
detected_electrode=0;
detected_electrode_predictive=0;

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
    for number=1:length(crisis{p})
        onset_electrodes=electrode_onset{p};
        binary_signal_2{crisis{p}(number)}=zeros(length(onset_electrodes{number}),length(init_register(number):init_register(number+1)));
        binary_signal{crisis{p}(number)}=zeros(1,length(init_register(number):init_register(number+1)));
        binary_signal{crisis{p}(number)}(1,electrodo_init(number):electrodo_end(number))=1;
        for j=1:length(onset_electrodes{number})
             sort_value=sort(test_distribution(onset_electrodes{number}(j),:));
             test_value=sort_value(900);
             index=find(all_eigenvalues(onset_electrodes{number}(j),floor(init_register(number)):floor(init_register(number+1)))>test_value);
             binary_signal_2{crisis{p}(number)}(j,index)=1;
        end
    end
    
    for k=1:length(binary_signal)
        seizure_flag=0;
        predicted_flag=0;
        for j=1:size(binary_signal_2{k},1)
            if isempty(binary_signal_2{k})~=1 && isempty(binary_signal{k})~=1       
                s_real=find(binary_signal{k});
                s_detected=find(binary_signal_2{k}(j,:));
                equal_data=0;
                for l=1:length(s_real)
                    for m=1:length(s_detected)
                        if s_real(l)==s_detected(m)+50
                            equal_data=equal_data+1;
                        end
                    end
                end
                if equal_data
                    detected_electrode=detected_electrode+1;
                    seizure_flag=1;
                else
                    if max(s_detected)+50<min(s_real)
                        detected_electrode_predictive=detected_electrode_predictive+1;
                        predicted_flag=1;
                    end
                end
                counter=counter+1;
            end 
        end
        detected_seizure=detected_seizure+seizure_flag;
        predicted_seizure=predicted_seizure+predicted_flag;
    end
    clear binary_signal binary_signal_2 electrodo_init electrodo_end init_register last_end
end

performance_electrode_detection=detected_electrode/counter;
performance_seizure_detection=detected_seizure/total_seizure;