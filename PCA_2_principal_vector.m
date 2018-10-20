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
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/PCA2/Histograms/';

crisis={[2,3,5,6,7],[2,3],[1,3,4]};
bands=['delta';'theta';'alpha';'beta_';'gamma'];

for p=1:3
    for number=1:length(crisis{p})
        patient_name=patients{p}.name;
        day=crisis{p}(number);
        name=strcat(patient_name,'_',int2str(day));
        load(name,'-mat')
        c = colormap(jet(length(find(pvalue{p,number}))));
        counter=1;
        for i=1:length(vector_projection{p}(number,:,:,:))-2
            significant_pos=find(significance{p}(number,i,:));
            if length(significant_pos)==1 && ~isempty(find(significant_pos==5))
                vector=reshape(vector_projection{p}(number,i,:,significant_pos),[5 1]);
                vector=vector/norm(vector);
                plot(1:5,abs(vector),'-o','Color',c(counter,:))
                hold on
                channels(counter,:)=header.channels(i,5:8);
                counter=counter+1;
            end 
        end
        xticks([1:5])
        xticklabels({bands})
        legend({channels})
        xtickangle(90)
        xlabel('Band Frequency')
        ylabel('Projection')
        box on
    end
end


