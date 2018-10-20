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
        new_index=x1(length(x1)-round(length(x1)*percent/100):end);
        figure(1)
        c = colormap(jet(length(new_index)));
        for indice=1:length(new_index)
                i=new_index(indice);
                vector=reshape(vector_projection{p}(number,i,:,5),[5 1]);
                vector=vector/norm(vector);
                plot(1:5,abs(vector),'-o','Color',c(indice,:))
                hold on
                channels(indice,:)=header.channels(i,5:8);
        end
        xticks([1:5])
        xticklabels({bands})
        legend({channels})
        xtickangle(90)
        xlabel('Band Frequency')
        ylabel('Projection')
        box on
        name=strcat(directory,patients{p}.name,'_',int2str(crisis{p}(number)),'_PCA2');
        saveas(1,name,'png')
        close(1)
    end
end
