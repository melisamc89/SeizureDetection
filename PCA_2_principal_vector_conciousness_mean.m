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

percent=30;
counter1=1;
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
        [val index]=sort(pvalue{p,number});
        x1=index(find(val));
        new_index=x1(length(x1)-round(length(x1)*percent/100)+1:end);
        counter=1;
        for indice=1:length(new_index)
                i=new_index(indice);
                vector=reshape(vector_projection{p}(number,i,:,5),[5 1]);
                vector=vector/norm(vector);
                band_analysis(counter,:)=vector;
                counter=counter+1;
        end
        band_analysis_mean{p}(number,:)=mean(abs(band_analysis));
        band_analysis_std{p}(number,:)=std(abs(band_analysis));
        counter1=counter1+1;
        clear band_analysis
    end
end

for p=1:3
    counter1=1;
    for number=1:length(crisis{p})
        loss{p}(counter1)=patients{p}.CSS(crisis{p}(number));
        counter1=counter1+1;
    end
end

bands=['delta';'theta';'alpha';'beta_';'gamma'];

for p=1:3
    for i=1:5
        subplot(2,3,i)
        errorbar(loss{p},band_analysis_mean{p}(:,i),band_analysis_std{p}(:,i),'o')
        hold on
        scatter(loss{p},band_analysis_mean{p}(:,i),'filled')
        box on
        xlabel('CSS')
        ylabel('proj')
        title(bands(i,:))
    end
end

%% for regresion plot
counter=1;
for p=1:3
    for i=1:length(crisis{p})
        bands_analysis_new(counter,:)=band_analysis_mean{p}(i,:);
        loss_new(counter)=loss{p}(i);
        counter=counter+1;
    end
end

for j=1:5
    aux=corrcoef(bands_analysis_new(:,j),loss_new);
    corr_coeff(j)=aux(1,2);
    [r m b]=regression(loss_new,bands_analysis_new(:,j)');
    r_val(j)=r;
    m_val(j)=m;
    b_val(j)=b;
end


bands=['delta';'theta';'alpha';'beta_';'gamma'];

for p=1:3
    for i=1:5
        subplot(2,3,i)
        errorbar(loss{p},band_analysis_mean{p}(:,i),band_analysis_std{p}(:,i),'o')
        hold on
        scatter(loss{p},band_analysis_mean{p}(:,i),'filled')
        box on
        xlabel('CSS')
        ylabel('proj')
        title(strcat(bands(i,:),':',num2str(corr_coeff(i))))
        
    end
end
for i=1:5
    x=0:0.1:10;
    line=m_val(i)*x+b_val(i);
    subplot(2,3,i)
    plot(x,line)
end
