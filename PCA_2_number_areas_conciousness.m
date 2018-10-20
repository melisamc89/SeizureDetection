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
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/PCA2/Eigenvectors/';

crisis={[2,3,5,6,7],[2,3],[1,3,4]};
bands=['delta';'theta';'alpha';'beta_';'gamma'];


pval_limit=0.01;
p_values=[0.01,0.001,0.0001,0.00001,0.000001,0.0000001];
for j=1:6
    pval_limit=p_values(j);
    counter2=1;
    for p=1:3
        counter3=1;
        for number=1:length(crisis{p})
            patient_name=patients{p}.name;
            day=crisis{p}(number);
            name=strcat(patient_name,'_',int2str(day));
            load(name,'-mat')
            counter1=1;
            for i=1:length(vector_projection{p}(number,:,:,:))
                significant_pos=find(significance{p}(number,i,:));
                if length(significant_pos)==1 && ~isempty(find(significant_pos==5)) && pvalue{p,number}(i)<pval_limit
                        electrode_name(counter1,:)=header.channels(i,5:8);
                        counter1=counter1+1;
                end
            end
            [involved_areas number_electrodes_per_area]=Involved_Areas(electrode_name);        
            loss(counter2)=patients{p}.CSS(crisis{p}(number));
            area_number(counter2)=length(involved_areas);
            counter2=counter2+1;
            loss_patient{p}(counter3)=patients{p}.CSS(crisis{p}(number));
            number_electrodes_area{p}(counter3)=length(involved_areas);
            counter3=counter3+1;
            clear electrode_name
          end
    end

    figure(j)
    aux=corrcoef(loss,area_number);
    corr_coeff=aux(1,2);
    [r m b]=regression(loss,area_number);

    for p=1:3
        scatter(loss_patient{p},number_electrodes_area{p},'filled');
        hold on
    end
    x=0:0.1:10;
    line=m*x+b;
    plot(x,line)
    box on
    legend('wagner','rizzio','molina','fit')
    xlabel('CSS')
    ylabel(strcat('Number of areas p_v_a_l < ', num2str(pval_limit)))
    title(strcat('CorrCoeff:',num2str(corr_coeff)))
    file_name=strcat(directory,'number_areas_vs_CSS_regression_p=',num2str(pval_limit));
    saveas(j,file_name,'png')
    close(j)
end


%% Idem but for CorrCoeff vs pvalue_limit

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

p_values=[0.000001,0.00001,0.0001,0.001,0.01,0.1,1];
for j=1:length(p_values)
    pval_limit=p_values(j);
    counter2=1; 
    for p=1:3
        counter3=1;
        for number=1:length(crisis{p})
            patient_name=patients{p}.name;
            day=crisis{p}(number);
            name=strcat(patient_name,'_',int2str(day));
            load(name,'-mat')
            counter1=1;
            electrode_name=[];
                for i=1:length(vector_projection{p}(number,:,:,:))
                    significant_pos=find(significance{p}(number,i,:));
                    if length(significant_pos)==1 && ~isempty(find(significant_pos==5)) && pvalue{p,number}(i)<pval_limit
                            electrode_name(counter1,:)=header.channels(i,5:8);
                            counter1=counter1+1;
                    end
                end
                if length(electrode_name)
                    [involved_areas number_electrodes_per_area]=Involved_Areas(electrode_name);        
                    loss(counter2)=patients{p}.CSS(crisis{p}(number));
                    area_number(counter2)=length(involved_areas);
                    counter2=counter2+1;
                    loss_patient{p}(counter3)=patients{p}.CSS(crisis{p}(number));
                    number_electrodes_area{p}(counter3)=length(involved_areas);
                    counter3=counter3+1;
                    clear electrode_name
                end
          end
    end 
    aux=corrcoef(loss,area_number);
    corr_coeff(j)=aux(1,2);
end

scatter(p_values(1:length(corr_coeff)),corr_coeff,'filled')
grid on
box on
xlabel('P valor limit')
ylabel('PCC')
set(gca, 'XScale', 'log')

