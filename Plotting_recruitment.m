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
electrode_onset_physician={{[1,2,3,13,14,15],[4,5,6],[13,14],[14],[3],[2,3,13,14],[3,4],[11,12,13],[11,12,13]},...
    {[25,26,26],[25,26,27],[25,26],[25,26],[26],[25,26]},...
    {[10,11,19,20],[28,29,37,38],[28,29,37,38,46,47],[10,11,19,20]},...
    {[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6]},...
    {[1,2,3,19,10,11],[37,46,47],[28,29,37,38,46,47],[19,20,10,11],[19,10,11]}};

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
            mark2=mark1(find(diff(mark1)>1));   
            if isempty(mark1)
                 index(i,number)=0;
            else
                if isempty(mark2)
                    mark=min(mark1);
                    index(i,number)=mark+floor(init_register(number))+1;
                else
                    auxiliar=mark2;
                    mark2(1)=mark1(1);
                    mark2(2:length(auxiliar)+1)=auxiliar;
                    mark2(end+1)=mark1(end);
                    sumador=zeros(1,length(mark2)-1);
                    for s=1:length(mark1)
                        for t=1:length(mark2)-1
                            if mark1(s)>mark2(t) && mark1(s)<mark2(t+1)
                                sumador(t)=sumador(t)+1;
                            end
                        end
                    end
                    mark3=mark2(find(sumador>51));
                    mark=max(mark3);
                    if isempty(mark)
                        index(i,number)=0;
                    else
                        index(i,number)=mark+floor(init_register(number))+1;
                    end
                end
            end
        end
    end
     
    for number=1:length(crisis{p})
        [value_1 order_1]=sort(index(:,number));
        value=value_1(find(value_1));
        order=order_1(find(value_1));
        init_channels=electrode_onset_physician{p}{number};
        init_channel_order=[];
        counter_init=1;
        counter_init_2=1;
        n=length(order);
        for j=1:length(init_channels)
            x=find(order==init_channels(j));
            if x
            init_channel_order(counter_init)=x;
            counter_init=counter_init+1;
            else
            order(n+counter_init_2)=init_channels(j);
            value(n+counter_init_2)=all_electrodo_end(number);
            init_channel_order(counter_init)=n+counter_init_2;
            counter_init=counter_init+1;
            counter_init_2=counter_init_2+1;
            end
        end
        
        %if isempty(init_channel_order) || length(init_channels_order)<length(init_channels)   
           % order=[order;init_channels'];
           % value=[value;ones(length(init_channels),1)*all_electrodo_end(number)'];
           % for j=1:length(init_channels)
           % x=find(order==init_channels(j));
           % init_channel_order(j)=x;
           % end
        %end
        order_channels=channels(order,:);
        inv_areas=Involved_Areas(channels(:,1:4));
        
        for l=1:length(inv_areas)
             for k=1:size(order_channels,1)
                 if strcmp(order_channels(k,1:3),inv_areas(l,:))
                    new_value(l,k)=value(k)+crisis_time;
                 end
             end
        end
        c=colormap(jet((length(inv_areas))*10));
        figure(1)
        %bar(1:length(w),1./w,'Barwidth',0.8)
        hold on
        for l=1:size(new_value,1)
             %scatter(1:length(header.channels),new_value(l,:),[],c((l-1)*10+1,:,:),'Filled')
             scatter(1:length(new_value(l,:)),new_value(l,:),[],c((l-1)*10+1,:,:),'Filled')
        end
        hold on
        scatter(init_channel_order,ones(length(init_channel_order),1)*all_electrodo_init(number),[],[0 0 0],'Filled')
        legend_areas=[inv_areas(1:end,:);'PhM'];
        legend({legend_areas})
        box on
        hold on
        line([0 length(order_channels)],[all_electrodo_init(number) all_electrodo_init(number)],'color',[1 0 0],'Linewidth',2)
        line([0 length(order_channels)],[all_electrodo_end(number) all_electrodo_end(number)],'color',[0 0 1],'Linewidth',2)
        xlabel('Channels')
        ylabel('Time [s]')
        xticks([1:length(order_channels)])
        xtickangle(90)
        xticklabels({order_channels})
        limit1=value(find(value));
        limit1(end+1)=all_electrodo_init(number);
        min_value=min(limit1)-5;
        limit2=value;
        limit2(end+1)=all_electrodo_end(number);
        max_value=max(limit2)+1;
        if isempty(min_value) min_value=-1; max_value=1; end
        ylim([min_value max_value+50])
        xlim([0 length(order_channels)])
        name=strcat(directory,patients(p,:),'_',num2str(crisis{p}(number)),'_Recluitment')
        saveas(1,name,'png')
        close(1)
        clear new_value order_channels
    end
    clear index
end

%% plot eigenvector
c=colormap(jet(50));
bands_projection=['D','T','A','B','G'];
for i=1:length(header.channels)
    figure(i)
    subplot(6,1,1)
    plot(all_eigenvalues(i,:),'b','Linewidth',2)
    hold on
    baseline=max(test_distribution(i,:));
    plot(baseline*ones(1,length(all_eigenvalues)),'r','Linewidth',2)
    box on
    ylabel('Eigenvalue')
    hold on
    for number=1:length(crisis{p})
        line([init_register(number) init_register(number)],[min(all_eigenvalues(i,:)) max(all_eigenvalues(i,:))],'Color','k')
        %hold on
        %line([electrodo_end(number)/10 electrodo_end(number)/10],[min(eigenvalue) max(eigenvalue)],'Color','b')
    end
    for l=1:5
        subplot(6,1,l+1)
        plot(abs(reshape(V_crisis_prueba(i,l,:),[length(V_crisis_prueba(i,l,:)),1])),'Color',c((l-1)*10+1,:,:),'Linewidth',2)
        box on
        %xlabel('Time [s]')
        ylabel(strcat('Proj:',bands_projection(l)))
        hold on
        for number=1:length(crisis{p})
            line([init_register(number) init_register(number)],[0 1],'Color','k')
        end
    end
    subplot(6,1,6)
    xlabel('Time [s]')
    nombre=strcat(directory,patient_name,'_',header.channels(i,5:8),'_eigenvector');
    saveas(i,nombre,'png')
    close(i)
end

%% plot eigenvalue

for i=1:length(header.channels)
    figure(i)
    plot(all_eigenvalues(i,:),'b','Linewidth',2)
    hold on
    baseline=max(test_distribution(i,:));
    plot(baseline*ones(1,length(all_eigenvalues)),'r','Linewidth',2)
    box on
    xlabel('Time [s]')
    ylabel('Eigenvalue')
    hold on
    for number=1:length(crisis{p})
        line([init_register(number) init_register(number)],[min(all_eigenvalues(i,:)) max(all_eigenvalues(i,:))],'Color','k')
        %hold on
        %line([electrodo_end(number)/10 electrodo_end(number)/10],[min(eigenvalue) max(eigenvalue)],'Color','b')
    end
    title(header.channels(i,5:8))
    nombre=strcat(directory,patient_name,'_',header.channels(i,5:8));
    saveas(i,nombre,'png')
    close(i)
end

%% plot signal projection


for i=1:length(header.channels)-2
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
        data2=(D^-0.5)*V'*(matrix-repmat(mean(matrix,2),1,b));     
        
        prueba=dot(reshape(V_crisis_prueba(i,:,:),[size(V_crisis_prueba(i,:,:),2) size(V_crisis_prueba(i,:,:),3)]),...
            (matrix(:,1:size(V_crisis_prueba(i,:,:),3))-repmat(mean(matrix,2),1,size(V_crisis_prueba(i,:,:),3))));                
        figure(i)
        subplot(2,1,1)
        plot(crisis_time:length(all_eigenvalues)+crisis_time-1,all_eigenvalues(i,:),'b','Linewidth',2)
        hold on
        title(header.channels(i,5:8))
        baseline=max(test_distribution(i,:));
        plot(baseline*ones(1,length(all_eigenvalues)),'r','Linewidth',2)
        box on
        ylabel('Eigenvalue')
        hold on
        for number=1:length(crisis{p})
            line([init_register(number) init_register(number)],[min(all_eigenvalues(i,:)) max(all_eigenvalues(i,:))],'Color','k')
            %hold on
            %line([electrodo_end(number)/10 electrodo_end(number)/10],[min(eigenvalue) max(eigenvalue)],'Color','b')
        end
        subplot(2,1,2)
        plot(abs(prueba),'k','Linewidth',2)
        hold on
        for number=1:length(crisis{p})
            line([init_register(number) init_register(number)],[min(abs(prueba)) max(abs(prueba))],'Color','r')
            %hold on
            %line([electrodo_end(number)/10 electrodo_end(number)/10],[min(eigenvalue) max(eigenvalue)],'Color','b')
        end
        box on
        ylabel('PCA')
        xlabel('Time [s]')
        nombre=strcat(directory,patient_name,'_',header.channels(i,5:8),'_projection');
        saveas(i,nombre,'png')
        close(i)
end

%% index 

for i=1:length(header.channels)-2
        baseline=max(test_distribution(i,:));
        for number=1:length(crisis{p})
            recording=all_eigenvalues(i,floor(init_register(number)):floor(init_register(number+1)));
            mark=find(recording>=baseline,1);
            if isempty(mark)
                index(i,number)=0;
            else 
                index(i,number)=mark+floor(init_register(number));
            end
        end
end
