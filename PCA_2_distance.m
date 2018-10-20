clear all
close all
data_base
load('eig_significance2','-mat');
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/PCA2/';
Time=1;
crisis={[2,3,5,6,7],[2,3],[1,3,4]};

%for electrodes
for p=1:3
    for number=1:length(crisis{p})
        patient_name=patients{p}.name;
        day=crisis{p}(number);
        count=1;
        for i=1:length(vector_projection{p}(number,:,:,:))
            significant_pos=find(significance{p}(number,i,:));
            if length(significant_pos)==1 && ~isempty(find(significant_pos==5))
                significant_vector_for_distance(:,count)=reshape(vector_projection{p}(number,i,:,significant_pos),[5 1])/...
                    norm(reshape(vector_projection{p}(number,i,:,significant_pos),[5 1]));
                count=count+1;
            end
        end
        count2=1;
        for j=1:count-1
            for k=j+1:count-1
                x=dot(significant_vector_for_distance(:,j),significant_vector_for_distance(:,k))/...
                    (norm(significant_vector_for_distance(:,j))*norm(significant_vector_for_distance(:,k)));
                dist{p,number}(count2)=atan(sqrt(1-x*x)/x);
                count2=count2+1;
            end
        end
        clear significant_vector_for_distance
    end
    %clear significant_vector_for_distance
end

subplot(1,2,1)
for i=1:length(significant_vector_for_distance)
   plot3([0 significant_vector_for_distance(1,i)],[0 significant_vector_for_distance(2,i)],[0 significant_vector_for_distance(3,i)])
   hold on
end
subplot(1,2,2)
for i=1:length(significant_vector_for_distance)
   plot3([0 significant_vector_for_distance(1,i)],[0 significant_vector_for_distance(4,i)],[0 significant_vector_for_distance(5,i)])
   hold on
end

for p=1:3
    for number=1:length(crisis{p})
        a=histogram(dist{p,number},25,'Normalization','probability')
        hold on
    end
end

%% for crisis

%for electrodes
for p=1:3
    electrode=zeros(length(crisis{p}),60);
    for number=1:length(crisis{p})
        patient_name=patients{p}.name;
        day=crisis{p}(number);
        for i=1:length(vector_projection{p}(number,:,:,:))
            significant_pos=find(significance{p}(number,i,:));
            if length(significant_pos)==1 && ~isempty(find(significant_pos==5))
                electrode(number,i)=1;
            end
        end
    end
    sig_electrodes=find(sum(electrode)==length(crisis{p}));
    for number=1:length(crisis{p})
        patient_name=patients{p}.name;
        day=crisis{p}(number);
        for i=1:length(sig_electrodes)
            significant_vector_for_distance(number,i,:)=reshape(vector_projection{p}(number,sig_electrodes(i),:,5),[5 1]);
        end
    end
    
    for j=1:length(sig_electrodes)
        count=1;
        for k=1:length(crisis{p})
            for l=k+1:length(crisis{p})
                x1=reshape(significant_vector_for_distance(k,j,:),[length(significant_vector_for_distance(k,j,:)) 1]);
                x2=reshape(significant_vector_for_distance(l,j,:),[length(significant_vector_for_distance(l,j,:)) 1]);
                x=dot(x1,x2)/...
                    (norm(x1)*norm(x2));
                dist{p}(j,count)=atan(sqrt(1-x*x)/x);
                count=count+1;
            end
        end
    end
    
    clear electrode significant_vector_for_distance
end
