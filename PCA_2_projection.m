clear all
close all
directory1='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/';
cd(directory1)
data_base
directory2='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/PCA2/';
cd(directory2)
load('eig_significance2','-mat');
cd(directory1)
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/PCA2/Histograms/Bin=500/';
Frequencies=0:0.1:100;
bands=[1,3.5,7.4,12.4,24,97];
Time=1;

crisis={[2,3,5,6,7],[2,3],[1,3,4]};
bin_width=500;

for p=1:3
    for number=1:length(crisis{p})
        patient_name=patients{p}.name;
        day=crisis{p}(number);
        name=strcat(patient_name,'_',int2str(day));
        load(name,'-mat')
        s_rate=header.samplingrate;
        Fs = s_rate;
        electrodo_init=patients{p}.crisis.electrodo_init(day)/s_rate;
        electrodo_end=patients{p}.crisis.electrodo_end(day)/s_rate;
        hipotesis{p,number}=zeros(length(vector_projection{p}(number,:,:,:)),1);
        pvalue{p,number}=zeros(length(vector_projection{p}(number,:,:,:)),1);
        for i=1:length(vector_projection{p}(number,:,:,:))
            %significant_pos=find(significance{p}(number,i,:));
            significant_pos=5;
            if length(significant_pos)==1 && ~isempty(find(significant_pos==5))
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
                %significant_projection(count,:,:)=reshape(vector_projection{p}...
                %(number,i,:,significant_pos(find(significant_pos==5))),...
                    %[5 1])'*(matrix-repmat(mean(matrix,2),1,b))/...
                   % norm(reshape(vector_projection{p}(number,i,:,significant_pos(find(significant_pos==5))),[5 1]));
                    prueba=reshape(vector_projection{p}(number,i,:,significant_pos(find(significant_pos==5))),...
                    [5 1])'*(matrix-repmat(mean(matrix,2),1,b))/...
                    norm(reshape(vector_projection{p}(number,i,:,significant_pos(find(significant_pos==5))),[5 1]));
                %plot(prueba)
                %hold on
                %count=count+1;
                tiempo=0:Time/2:length(matrix)*Time/2-Time/2;
                count1=1;
                count2=1;
                for j=1:length(tiempo)
                    if tiempo(j)>electrodo_init && tiempo(j)<electrodo_end
                        data1(count1)=prueba(j);
                        count1=count1+1;
                    else
                        data2(count2)=prueba(j);
                        count2=count2+1;
                    end
                end
                %figure(i)
                %a = histogram(data1,'BinWidth',bin_width,'Normalization','probability');
                %hold on
                %b = histogram(data2,'BinWidth',bin_width,'Normalization','probability');
                %hold off
                %alpha(b,.5)
                %xlabel('Projection Value')
                %ylabel('Probability')
                %legend('Ictal','InterIctal')
                %title(strcat(patient_name,':',int2str(day),':',header.channels(i,5:8)))
                %name=strcat(directory,patient_name,'_',int2str(day),'_',header.channels(i,5:8));
                %saveas(i,name,'png')
                %close(i)
                figure(i)
                a = histogram(data1,'BinWidth',bin_width,'Normalization','cdf');
                hold on
                b = histogram(data2,'BinWidth',bin_width,'Normalization','cdf');
                %xlabel('Projection Value')
                %ylabel('Probability')
                %legend('Ictal','InterIctal')
                %title(strcat(patient_name,':',int2str(day),':',header.channels(i,5:8)))
                %name=strcat(directory,patient_name,'_',int2str(day),'_',header.channels(i,5:8),'_cdf');
                %saveas(i,name,'png')
                aux1=a.Values;%(1:min(length(a.Values),length(b.Values)));
                aux2=b.Values;%(1:min(length(a.Values),length(b.Values)));
                [h1,p1,kv1] = kstest2(aux1,aux2,'Alpha',0.01);
                hipotesis{p,number}(i)=h1;
                pvalue{p,number}(i)=p1;
                if p1==0
                    pvalue{p,number}(i)=eps;
                end
                close(i)
                clear data1 data2
            end 
        end
    end
end

directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/PCA2/';
nombre=strcat(directory,'hipotesis')
save(nombre,'hipotesis','pvalue','-mat')


