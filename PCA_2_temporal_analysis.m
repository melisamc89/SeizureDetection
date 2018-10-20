clear all
close all
directory1='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/';
cd(directory1)
data_base
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/PCA2/';
cd(directory)
load('hipotesis','-mat')
load('eig_significance2','-mat');
cd(directory1)
crisis={[2,3,5,6,7],[2,3],[1,3,4]};
Frequencies=0:0.1:100;
bands=[1,3.5,7.4,12.4,24,97];
Time=1;
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/PCA2/Pvalue/';

percent=20;
n=50;
bin_width=1000;

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
        electrodo_prop=patients{p}.crisis.propagation(day)/s_rate-electrodo_init;
        if electrodo_prop <0 electrodo_prop=0; end
        electrodo_AOC=patients{p}.crisis.AOC(day)/s_rate-electrodo_init;
        if electrodo_AOC < 0 electrodo_AOC=0; end
        [val index]=sort(pvalue{p,number});
        x1=index(find(val));
        new_index=index(length(index)-round(length(index)*percent/100):end);
        for indice=1:length(new_index)
            i=new_index(indice);
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
            significant_pos=5;
            prueba=reshape(vector_projection{p}(number,i,:,significant_pos(find(significant_pos==5))),...
                        [5 1])'*(matrix-repmat(mean(matrix,2),1,b))/...
                        norm(reshape(vector_projection{p}(number,i,:,significant_pos(find(significant_pos==5))),[5 1]));
           tiempo=0:Time/2:length(matrix)*Time/2-Time/2;            
           for k=1:n*2
                init=electrodo_init-(n-k)*Time/2;
                if init<0 init=0; end
                count1=1;
                count2=1;
                for j=1:length(tiempo)
                    if tiempo(j)>init && tiempo(j)<electrodo_end
                        data1(count1)=prueba(j);
                        count1=count1+1;
                    end
                    if tiempo(j)<init || tiempo(j)>electrodo_end;
                        data2(count2)=prueba(j);
                        count2=count2+1;
                    end
                end
                figure(1)
                a = histogram(data1,'BinWidth',bin_width,'Normalization','cdf');
                hold on
                b = histogram(data2,'BinWidth',bin_width,'Normalization','cdf');
                aux1=a.Values;%
                aux2=b.Values;%(1:min(length(a.Values),length(b.Values)));
                [h1,p1,kv1] = kstest2(aux1,aux2,'Alpha',0.01);
                close(1)
                hipotesis_varios(indice,k)=h1;
                pvalue_varios(indice,k)=p1;
                clear data1 data2
            end   
        end
        figure(1)
        channels=header.channels(:,5:8);
        used_channels=channels(new_index,:);
        counter=1;
        for i=1:10
        y=find(pvalue_varios(i,:)<0.1,1);
            if isempty(y)==0
            x(i)=y;
            else
            x(i)=1;
            end
        end
        numbers=[-n+1:n]/2;
        z=numbers(x);
        [value indice]=sort(z);
        used_channels(indice,:);
        subplot(2,1,1)
        c = colormap(jet(length(value)));        
        for k=1:length(value)
            plot([-n+1:n]/2,pvalue_varios(indice(k),:),'Color',c(k,:),'Linewidth',1.5)
            hold on
        end
        line([0 0],[0 max(max(pvalue_varios))],'Color','k','Linewidth',2)
        hold on
        line([electrodo_prop electrodo_prop],[0 max(max(pvalue_varios))],'Color','b','Linewidth',2)
        line([electrodo_AOC electrodo_AOC],[0 max(max(pvalue_varios))],'Color','r','Linewidth',2)
        legend({used_channels(indice,:)})
        xlabel('Time from Nuria s mark [s]')
        ylabel('pvalue')
        subplot(2,1,2)
        line([0 0],[0 max(max(pvalue_varios))],'Color','k','Linewidth',2)
        hold on
        line([electrodo_prop electrodo_prop],[0 max(max(pvalue_varios))],'Color','b','Linewidth',2)
        line([electrodo_AOC electrodo_AOC],[0 max(max(pvalue_varios))],'Color','r','Linewidth',2)
        title('Electrodes with pval<0.1')
        for k=1:length(value)
            line([value(k) value(k)],[0 max(max(pvalue_varios))],'Color',c(k,:),'Linewidth',1)
            hold on
        end
        box on
        %legend({used_channels(indice,:)})
        xlabel('Time from Nuria s mark [s]')
        ylabel('-')
        name=strcat(directory,patient_name,'_',int2str(crisis{p}(number)),'_pvalues_evolution2');
        saveas(1,name,'png')
        close(1)
        clear hipotesis_varios pvalue_varios
    end
end





