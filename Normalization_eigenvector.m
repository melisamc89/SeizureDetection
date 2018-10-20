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
sampling_rate=[2048,2000,2000,2000,2000];
crisis_time=25;

for p=1:5
    last_end=0;
    signal=[];
    for number=1:length(crisis{p})
        patient_name=patients(p,:);
        day=crisis{p}(number);
        name=strcat(directory2,patient_name,'_',int2str(day));
        load(name,'-mat')
        if p==1
            data1=data(1:44,:);clear data;data=data1;end
        s_rate=sampling_rate(p);
        electrodo_init(number)=SzO/s_rate+last_end;
        electrodo_end(number)=SzE/s_rate+last_end;
        init_register(number)=last_end;
        last_end=last_end+length(data(1,:))/s_rate;
        signal=[signal,data];
    end
    for i=1:size(signal,1)
        auxiliar=signal(i,:);
        window_size = s_rate*Time;
        overlap = 0;
        [Y,F,T,PSD] = spectrogram(auxiliar,window_size,overlap,Frequencies,s_rate,'power','yaxis');
        delta=sum(PSD(1*10:3.5*10,:));
        theta=sum(PSD(3.5*10:7.4*10,:));
        alpha=sum(PSD(7.4*10:12.4*10,:));
        beta=sum(PSD(12.4*10:24*10,:));
        gamma=sum(PSD(24*10:97*10,:));
        matrix=[delta;theta;alpha;beta;gamma];
        [a b]=size(matrix);    
        matrix2=(matrix-repmat(mean(matrix,2),1,b));
        tiempo=0:Time:length(matrix2)*Time-Time;            
        [V,D]=eig((matrix2*matrix2')/(length(matrix2)-1));
        max_norm_eigenvector(i,:,:)=V;
    end
    for i=1:5
        subplot(3,2,i)
        plot(abs(reshape(max_norm_eigenvector(:,:,i),[44 5])'))
    end
end


