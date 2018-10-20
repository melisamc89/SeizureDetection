clear all
close all
directory1='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/';
cd(directory1)
data_base
directory2='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/PCA2/';
cd(directory2)
load('eig_significance2','-mat');
cd(directory1)
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/PCA2/Pvalue/';
crisis={[2,3,5,6,7],[2,3],[1,3,4]};
Frequencies=0:0.1:100;
bands=[1,3.5,7.4,12.4,24,97];
Time=1;

for p=1:3
  count2=1;
    count1=1;
    for number=1:length(crisis{p})
        patient_name=patients{p}.name;
        day=crisis{p}(number);
        name=strcat(patient_name,'_',int2str(day));
        load(name,'-mat')
        s_rate=header.samplingrate;
        Fs = s_rate;
        electrodo_init=patients{p}.crisis.electrodo_init(day)/s_rate;
        electrodo_end=patients{p}.crisis.electrodo_end(day)/s_rate;
        for i=1:length(vector_projection{p}(number,:,:,:))
            significant_pos=find(significance{p}(number,i,:));
            if length(significant_pos) && ~isempty(find(significant_pos==5))
                    signal=data(i,:);
                    window_size = s_rate*Time;
                    overlap = window_size/2;
                    [Y,F,T,PSD] = spectrogram(signal,window_size,overlap,Frequencies,s_rate,'power','yaxis');
                    delta=sum(PSD(1*10:3.5*10,:));
                    theta=sum(PSD(3.5*10:7.4*10,:));
                    alpha=sum(PSD(7.4*10:12.4*10,:));
                    beta=sum(PSD(12.4*10:24*10,:));
                    gamma=sum(PSD(24*10:97*10,:));
                    matrix=[delta;theta;alpha;beta;gamma];
                    [a b]=size(matrix);
                    prueba=reshape(vector_projection{p}(number,i,:,significant_pos(find(significant_pos==5))),...
                        [5 1])'*(matrix-repmat(mean(matrix,2),1,b))/...
                        norm(reshape(vector_projection{p}(number,i,:,significant_pos(find(significant_pos==5))),[5 1]));
                    tiempo=0:Time/2:length(matrix)*Time/2-Time/2;            
                    for j=1:length(tiempo)
                        if tiempo(j)>electrodo_init && tiempo(j)<electrodo_end
                               data1{p}(count1)=prueba(j);
                                count1=count1+1;
                        end
                        if tiempo(j)<electrodo_init || tiempo(j)>electrodo_end;
                            data2{p}(count2)=prueba(j);
                            count2=count2+1;
                        end
                   end 
            end
        end
    end
end

bin_width=100;
figure(1)
for p=1:3
    subplot(3,1,p)
    a=histogram(data1{p},'BinWidth',bin_width,'Normalization','probability');
    hold on
    b=histogram(data2{p},'BinWidth',bin_width,'Normalization','probability');
    xlim([-10000 10000])
end
close(1)