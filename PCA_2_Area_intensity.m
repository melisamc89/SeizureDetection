clear all
close all
directory1='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/';
cd(directory1)
data_base
directory2='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/PCA2/';
cd(directory2)
load('eig_significance2','-mat');
load('hipotesis','-mat')
cd(directory1)
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/PCA2/Conciousness/';

crisis={[2,3,5,6,7],[2,3],[1,3,4]};

for p=1:3
    patient_name=patients{p}.name;
    day=crisis{p}(1);
    name=strcat(patient_name,'_',int2str(day));
    load(name,'-mat')
    [s t]=AreasIndex(header.channels(1:end-2,:));
    [inv_areas number]=Involved_Areas(header.channels(:,5:8));
    areas{p}=inv_areas;
    for number=1:length(crisis{p})
        p_val=pvalue{p,number};
        for j=1:length(s)
            x=p_val(s(j):t(j));
            value=min(x(find(x)));
            if value
            val{p}(number,j)=value;
            else
            val{p}(number,j)=0;
            end
        end
    end
end

for p=1:3
    counter1=1;
    for number=1:length(crisis{p})
        loss{p}(counter1)=patients{p}.CSS(crisis{p}(number));
        counter1=counter1+1;
    end
end

for p=1:3
    subplot(3,1,p)
    bar(1./val{p})
    legend(areas{p})
    xticklabels({loss{p}})
    set(gca, 'YScale', 'log')
    xlabel('CSS')
    ylabel('1/pvalue')
    box on
end