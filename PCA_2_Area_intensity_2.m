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
bands=['delta';'theta';'alpha';'beta_';'gamma'];


pval_limit=0.01;
p_values=[0.1,0.01,0.001,0.0001,0.00001,0.000001,0.0000001];
for j=1:7
    pval_limit=p_values(j);
    figure(j)
    for p=1:3
        counter3=1;
        patient_name=patients{p}.name;
        day=crisis{p}(1);
        name=strcat(patient_name,'_',int2str(day));
        load(name,'-mat')
        name_areas=NameAreas(header.channels);
        number_area=zeros(length(name_areas),length(crisis{p}));
        for number=1:length(crisis{p})
            counter1=1;
            for i=1:length(vector_projection{p}(number,:,:,:))
                significant_pos=find(significance{p}(number,i,:));
                if length(significant_pos)==1 && ~isempty(find(significant_pos==5)) && pvalue{p,number}(i)<pval_limit
                        electrode_name(counter1,:)=header.channels(i,5:8);
                        counter1=counter1+1;
                end
            end
            number_electrodes_per_area=Electrodes_per_area(electrode_name,header.channels);
            number_area(:,number)=number_electrodes_per_area;
            
            loss_patient{p}(counter3)=patients{p}.CSS(crisis{p}(number));
            counter3=counter3+1;
            clear electrode_name
        end
        subplot(3,1,p)
        [x index]=sort(loss_patient{p});
        number_area_2=number_area(:,index);
        bar(number_area_2')
        xlabel('CSS')
        ylabel('#')
        legend(name_areas)
        xticklabels({x})
    end

    file_name=strcat(directory,'number_electrodes_per_area_p=',num2str(pval_limit));
    saveas(j,file_name,'png')
    close(j)
end
