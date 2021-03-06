% Construct binary signals to compare automatic detection with physicians
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
crisis={[2,3,4,5,6,7,9,11,12],[2,3,4,5,6,7],[1,2,3,4],[37,38,48,50,58,62,66,68,68],[1,3,6,7,8]};
AOC_value={[0,3,5,4,6,1,0,1,1],[9,9,7,-1,0,-1],[9,9,7,9],[5,5,8,5,5,5,5,-1],[1,-1,6,6,6]};
electrode_onset_physician={{[1,2,3,13,14,15],[4,5,6],[13,14],[14],[3],[2,3,13,14],[3,4],[11,12,13],[11,12,13]},...
    {[25,26,26],[25,26,27],[25,26],[25,26],[26],[25,26]},...
    {[10,11,19,20],[28,29,37,38],[28,29,37,38,46,47],[10,11,19,20]},...
    {[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6],[3,4,5,6]},...
    {[1,2,3,19,10,11],[37,46,47],[28,29,37,38,46,47],[19,20,10,11],[19,10,11]}};

electrode_onset_physician2={{[1,5,6],[2],[15],[5],[1],[1,5],[1,2],[4,5],[4,5]},...
    {[9],[9],[9],[9],[9],[9]},...
    {[4,7],[10,13],[10,13,16],[4,7]},...
    {[1,2],[1,2],[1,2],[1,2],[1,2],[1,2],[1,2],[1,2],[1,2]},...
    {[1,4,7],[13,16],[10,13,16],[4,7],[4,7]}};

areas_delimitation={{[1,2,3;4,5,6;7,8,8;9,10,11;12,13,14;15,16,17;18,19,20;21,22,23;24,25,26;27,28,29;30,31,32;33,34,35;36,37,38;39,40,41;42,43,44]},...
    {reshape(1:1:45,[3 15])'},{reshape(1:1:54,[3 18])'},...
    {[1,2,3;4,5,6;7,7,7;8,9,10;11,12,13;14,14,14;21,21,21;15,16,17;18,19,20;22,23,23;24,25,25;26,27,28;29,30,31;32,32,32;33,34,35;36,37,38;39,39,39]},...
    {reshape(1:1:54,[3 18])'}};

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
    name_eigen=strcat(directory1,patient_name,'_all_eigenanalysis');
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
    for i=1:length(areas_delimitation{p}{1})
        new_channel_name(i,:)=channels(areas_delimitation{p}{1}(i,1),:);
        for j=1:length(areas_delimitation{p}{1}(i,:))
            sorted_test=sort(test_distribution(areas_delimitation{p}{1}(i,j),:));
            baseline=sorted_test(950);
            for number=1:length(crisis{p})
                recording=all_eigenvalues(i,floor(init_register(number))+1:floor(init_register(number+1)));
                mark1=find(recording>baseline);
                in_line=[electrodo_init(number)-50:1:electrodo_init(number)+200];
                [tf,loc]=ismember(in_line,mark1);
                if sum(tf)~=0
                    [val idx]=max(recording(in_line(find(tf))));
                    mark(number,j)=idx+floor(init_register(number))+1;
                else
                    mark(number,j)=0;
                end
            end
        end
        for number=1:length(crisis{p})
            if length(find(mark(number,:)))>2
                index(i,number)=min(mark(number,:));
            else
                index(i,number)=0;
            end
        end
    end
     
    for number=1:length(crisis{p})
        [value_1 order_1]=sort(index(:,number));
        value=value_1(find(value_1));
        order=order_1(find(value_1));
        init_channels=electrode_onset_physician2{p}{number};
        init_channel_order=[];
        counter_init=1;
        counter_init_2=1;
        n=length(order);
        for j=1:length(electrode_onset_physician2{p}{number})
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
        order_channels=new_channel_name(order,:);
        inv_areas=Involved_Areas(new_channel_name(:,1:4));
        
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