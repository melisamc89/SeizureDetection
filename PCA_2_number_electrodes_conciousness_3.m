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
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/PCA2/Recluting_test/';

percent=99;
n=15;
bin_width=1000;
NTIME=30;

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
        [val index]=sort(pvalue{p,number}(1:end-2));
        x1=index(find(val));
        new_index=index(length(index)-round(length(index)*percent/100)+1:end);
        for indice=32:length(new_index)
            i=new_index(indice);
            signal=data(i,:);
            window_size = s_rate*Time;
            overlap = 0;
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
           tiempo=0:Time:length(matrix)*Time-Time;            
           for k=1:n*2
                init=electrodo_init-(n-k)*Time;
                if init<0 init=0; end
                count1=1;
                count2=1;
                for j=1:length(tiempo)
                    if tiempo(j)>init && tiempo(j)<init+NTIME
                        data1(count1)=prueba(j);
                        count1=count1+1;
                    end
                    if tiempo(j)<init && tiempo(j)>electrodo_init-NTIME;
                        data2(count2)=prueba(j);
                        count2=count2+1;
                    end
                end
                if exist('data1') && exist('data2')
                figure(1)
                a = histogram(data1,'BinWidth',bin_width,'Normalization','cdf');
                hold on
                b = histogram(data2,'BinWidth',bin_width,'Normalization','cdf');
                aux1=a.Values;%
                aux2=b.Values;%(1:min(length(a.Values),length(b.Values)));
                [h1,p1,kv1] = kstest2(aux1,aux2,'Alpha',0.01);
                close(1)
                hipotesis_varios{p}(indice,k)=h1;
                pvalue_varios{p}(indice,k)=p1;
                pvalue_test{p}(indice,k,:)=Pvalue_Shuffling(data1,data2,bin_width,100);
                clear data1 data2
                else
                hipotesis_varios{p}(indice,k)=0;
                pvalue_varios{p}(indice,k)=1;
                pvalue_test{p}(indice,k,:)=ones(100,1);
                end
            end   
        end
    end
end

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
        [val index]=sort(pvalue{p,number}(1:end-2));
        x1=index(find(val));
        new_index=index(length(index)-round(length(index)*percent/100)+1:end);
        
        channels=header.channels(:,5:8);
        used_channels=channels(new_index,:);
        for i=1:length(new_index)
            for k=1:n*2
                auxiliar1=sort(1./pvalue_test{p}(i,k,:));
                plimit=auxiliar1(95);
                if 1/pvalue_varios{p}(i,k)>plimit
                   ptest(i,k)=1;
                else
                   ptest(i,k)=0;
                end
            end
        end
        channels_counter=1;
        for i=1:length(new_index)
            [val auxiliar1]=find(ptest(i,:));
            auxiliar2=abs(auxiliar1-n);
            y=auxiliar1(min(find(auxiliar2==min(auxiliar2))));
            if isempty(y)==0
                x(i)=y;
                x2(i)=pvalue_varios{p}(i,y);
            else
                x(i)=NTIME+1;
                x2(i)=1;
            end  
            if isempty(y)==0
                significant_channels(channels_counter,:)=used_channels(i,:);
                channels_counter=channels_counter+1;
            end
        end
        numbers=[-n+1:n+2];
        z=[numbers(x),0,electrodo_prop,electrodo_AOC];
        uninvolved_electrodes=find(x2==1);
        z(uninvolved_electrodes)=max(n,max(electrodo_prop,electrodo_AOC))+1;
        events=[used_channels;'Eve1';'Eve2';'Eve3'];
        w=[x2,min(x2),min(x2),min(x2)];
        [value indice]=sort(z);
        events=events(indice,:);
        w=w(indice);
        inv_areas=Involved_Areas(header.channels(:,5:8));
        inv_areas=[inv_areas;'Eve'];
        
        count=ones(length(inv_areas),1);
        for l=1:length(inv_areas)
            new_value2{l}=zeros(length(events),1);
            time_value2{l}=zeros(length(events),1);
            for k=1:length(events)
                if strcmp(events(k,1:3),inv_areas(l,:))
                    new_value2{l}(k)=w(k);
                    time_value2{l}(k)=value(k);
                end
            end
        end
        
        c=colormap(jet(length(inv_areas)));
        figure(1)
        subplot(2,1,1)
        %bar(1:length(w),1./w,'Barwidth',0.8)
        hold on
        for l=1:length(inv_areas)-1
            %if isempty(new_index2{l})==0
                bar(1:length(events),1./new_value2{l},'Barwidth',0.45,'Facecolor',c(l,:,:))
            %end
        end
        bar(1:length(events),1./new_value2{l+1},'Barwidth',0.1,'Facecolor',c(l+1,:,:))
        xticks([1:length(w)])
        xticklabels({events})
        xtickangle(90)
        xlabel('Events')
        ylabel('1/pvalue')     
        set(gca, 'YScale', 'log')
        ylim([0 1/min(x2)])
        box on
        subplot(2,1,2)
        hold on
        for l=1:length(inv_areas)-1
            %if isempty(new_index2{l})==0
                bar(1:length(events),time_value2{l},'Barwidth',0.45,'Facecolor',c(l,:,:))
            %end
        end
        hold on
        %places=find(new_value2{l+1});
        %eventos=zeros(length(events),1);
        %eventos(places)=max(value);
        bar(1:length(events),time_value2{l+1},'Barwidth',0.1,'Facecolor',c(l+1,:,:))
        xticks([1:length(w)])
        xticklabels({events})
        xtickangle(90)
        xlabel('Events')
        ylabel('Time')
        ylim([min(z) max(z)])
        box on
        
        name=strcat(directory,patient_name,'_',int2str(crisis{p}(number)),'_pvalues_rec');
        saveas(1,name,'png')
        close(1)
        pre{p}(number)=length(find(value<=0));
        inter{p}(number)=length(find(value>0 & value<max(electrodo_prop,electrodo_AOC))); 
        clear hipotesis_varios pvalue_varios new_index2 new_value2 time_value2 x x2
    end
end

counter=1;
for p=1:3
    for number=1:length(crisis{p})
        loss_patient{p}(number)=patients{p}.CSS(crisis{p}(number));
        loss(counter)=patients{p}.CSS(crisis{p}(number));
        pre_all(counter)=pre{p}(number);
        inter_all(counter)=inter{p}(number);
        counter=counter+1;
    end
end

for p=1:3
    subplot(3,1,1)
    scatter(loss_patient{p},pre{p},'filled')
    hold on
    box on
    subplot(3,1,2)
    scatter(loss_patient{p},inter{p},'filled')
    hold on
    box on
    subplot(3,1,3)
    scatter(loss_patient{p},pre{p}+inter{p},'filled')
    hold on
    box on
end
legend('wagner','rizzio','molina')

[r1 m1 b1]=regression(loss,pre_all);
[r2 m2 b2]=regression(loss,inter_all);
[r3 m3 b3]=regression(loss,pre_all+inter_all);
x=0:0.1:10;
line1=m1*x+b1;
subplot(3,1,1)
plot(x,line1)
ylabel('#Electrodes')
title(strcat('Corr=',num2str(r1)))
line2=m2*x+b2;
subplot(3,1,2)
plot(x,line2)
ylabel('#Electrodes')
title(strcat('Corr=',num2str(r2)))
line3=m3*x+b3;
subplot(3,1,3)
plot(x,line3)
ylabel('#Electrodes')
xlabel('CSS')
title(strcat('Corr=',num2str(r2)))
