% Construct binary signals to compare automatic detection with physicians
% detection
clear all
close all
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/26-07-2018/';
directory2='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/Datos/Datos_MATLAB/';
directory3='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/26-07-2018/25s_SlidingWindow/Eigenvalues_95/';
cd(directory)
Frequencies=0:0.1:100;
bands=[1,3.5,7.4,12.4,24,97];
Time=1;
patients=['WL';'RF';'MA';'CG';'CM'];
crisis={[2,3,4,5,6,7,9,11,12],[2,3,4,5,6,7],[1,2,3,4],[37,38,48,50,58,62,66,68,68],[1,3,6,7,8]};
electrode_onset_physician={{[1,2,3,13,14,15],[4,5,6],[13,14],[14],[3],[2,3,13,14],[3,4],[11,12,13],[11,12,13]},...
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
    N=length(electrode_onset_physician{p});
    init_register=floor(init_register);
    for i=1:N
        %test=sort(test_distribution(:));
        %value=test(95*size(all_eigenvalues,1)*10);
        for elec=1:length(electrode_onset_physician{p}{i})
            figure(1)
            test=sort(test_distribution(electrode_onset_physician{p}{i}(elec),:));
            value=test(950);
            plot(init_register(i):init_register(i+1),all_eigenvalues(electrode_onset_physician{p}{i}(elec),init_register(i):init_register(i+1)))
            hold on
            plot(init_register(i):init_register(i+1),ones(length(init_register(i):init_register(i+1)),1)*value,'r')
            line([all_electrodo_init(i) all_electrodo_init(i)],[0 45],'color',[1 0 0],'Linewidth',2)
            line([all_electrodo_end(i) all_electrodo_end(i)],[0 45],'color',[0 0 1],'Linewidth',2)
            xlabel('Time (s)')
            ylabel('Eigenvalue')
            legend('Eigenvalue','Test','Onset','End')
            name=strcat(directory3,patients(p,:),'_',num2str(crisis{p}(i)),'_',num2str(electrode_onset_physician{p}{i}(elec)));
            saveas(1,name,'png')
            close(1)
        end
    end
end