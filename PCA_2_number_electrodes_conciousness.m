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


p_values=[0.01,0.001,0.0001,0.00001,0.000001,0.0000001,0.00000001];
for j=1:7
    pval_limit=p_values(j);
    counter3=1;
    electrode_number=zeros(10,1);
    for p=1:3
        counter2=1;
        electrodes_number_patient{p}=zeros(length(crisis{p}),1);
        for number=1:length(crisis{p})
            for i=1:length(vector_projection{p}(number,:,:,:))
                significant_pos=find(significance{p}(number,i,:));
                if length(significant_pos)==1 && ~isempty(find(significant_pos==5)) && pvalue{p,number}(i)<pval_limit
                        electrode_number(counter3)=electrode_number(counter3)+1;
                        electrodes_number_patient{p}(counter2)=electrodes_number_patient{p}(counter2)+1;
                end
            end
            loss(counter3)=patients{p}.CSS(crisis{p}(number));
            loss_patient{p}(counter2)=patients{p}.CSS(crisis{p}(number));
            counter2=counter2+1;
            counter3=counter3+1;
        end
    end
   
    figure(j)
    aux=corrcoef(loss,electrode_number);
    corr_coeff=aux(1,2);
    [r m b]=regression(loss,electrode_number');

    for p=1:3
        scatter(loss_patient{p},electrodes_number_patient{p},'filled');
        hold on
    end
    x=0:0.1:10;
    line=m*x+b;
    plot(x,line)
    box on
    legend('wagner','rizzio','molina','fit')
    xlabel('CSS')
    ylabel(strcat('Number of electrodes p_v_a_l < ', num2str(pval_limit)))
    title(strcat('CorrCoeff:',num2str(corr_coeff)))
    file_name=strcat(directory,'number_electrodes_vs_CSS_regression_p=',num2str(pval_limit));
    saveas(j,file_name,'png')
    close(j)
end



%% Idem before but for Corr_Coeff vs P_limit

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


p_values=[0.00000001:0.00000001:0.0000001,0.0000001:0.0000001:0.000001,0.000001:0.000001:0.00001,...
    0.00001:0.00001:0.0001,0.0001:0.0001:0.001,0.001:0.001:0.01,0.01:0.01:0.1,0.1:0.1:1];
for j=1:length(p_values)
    pval_limit=p_values(j);
    counter3=1;
    electrode_number=zeros(10,1);
    for p=1:3
        counter2=1;
        electrodes_number_patient{p}=zeros(length(crisis{p}),1);
        for number=1:length(crisis{p})
            for i=1:length(vector_projection{p}(number,:,:,:))
                significant_pos=find(significance{p}(number,i,:));
                if length(significant_pos)==1 && ~isempty(find(significant_pos==5)) && pvalue{p,number}(i)<pval_limit
                        electrode_number(counter3)=electrode_number(counter3)+1;
                        electrodes_number_patient{p}(counter2)=electrodes_number_patient{p}(counter2)+1;
                end
            end
            loss(counter3)=patients{p}.CSS(crisis{p}(number));
            loss_patient{p}(counter2)=patients{p}.CSS(crisis{p}(number));
            counter2=counter2+1;
            counter3=counter3+1;
        end
    end
    aux=corrcoef(loss,electrode_number);
    corr_coeff(j)=aux(1,2);
end

scatter(p_values,corr_coeff,'filled')
grid on
box on
xlabel('P valor limit')
ylabel('PCC')
set(gca, 'XScale', 'log')




