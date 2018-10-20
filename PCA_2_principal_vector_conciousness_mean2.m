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

p_values=[0.00000001:0.00000001:0.0000001,0.0000001:0.0000001:0.000001,0.000001:0.000001:0.00001,...
    0.00001:0.00001:0.0001,0.0001:0.0001:0.001,0.001:0.001:0.01,0.01:0.01:0.1,0.1:0.1:1];
for j=1:length(p_values)
    pval_limit=p_values(j);
    counter1=1;
    for p=1:3
        for number=1:length(crisis{p})
            val=pvalue{p,number};
            x1=val(find(val));
            index=find(val);
            new_index=index(find(x1<pval_limit));
            band_analysis=[];
            counter=1;
            for indice=1:length(new_index)
                    i=new_index(indice);
                    vector=reshape(vector_projection{p}(number,i,:,5),[5 1]);
                    vector=vector/norm(vector);
                    band_analysis(counter,:)=vector;
                    counter=counter+1;
            end
            if counter>1
            band_analysis_mean{p}(number,:)=mean(abs(band_analysis));
            counter1=counter1+1;
            clear band_analysis
            end
        end
    end
    counter2=1;
    for p=1:3
        counter2=1;
        for number=1:length(crisis{p})
            loss{p}(counter2)=patients{p}.CSS(crisis{p}(number));
            counter2=counter2+1;
        end
    end
    if counter1>1
        counter=1;
        for p=1:3
            for i=1:length(crisis{p})
                bands_analysis_new(counter,:)=band_analysis_mean{p}(i,:);
                loss_new(counter)=loss{p}(i);
                counter=counter+1;
            end
        end
        for k=1:5
            aux=corrcoef(bands_analysis_new(:,k),loss_new);
            corr_coeff(j,k)=aux(1,2);
        end
    end
end

for j=1:5
    subplot(2,3,j)
    plot(p_values,corr_coeff(:,j),'Linewidth',2)
    box on
    xlabel('plimit')
    ylabel('PCC')
    title(bands(j,:))
    set(gca, 'XScale', 'log')
    ylim([-1 1])
end
