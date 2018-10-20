clear all
close all
data_base

directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/PCA2/';
Frequencies=0:0.1:100;
bands=[1,3.5,7.4,12.4,24,97];
Time=1;
cannonic=[1,1,1,1,1];
crisis={[2,3,5,6,7],[2,3],[1,3,4]};

for p=1:3
    for number=1:length(crisis{p})
        patient_name=patients{p}.name;
        day=crisis{p}(number);
        name=strcat(patient_name,'_',int2str(day));
        load(name,'-mat')
        s_rate=header.samplingrate;
        Fs = s_rate;            % Sampling frequency
        electrodo_init=patients{p}.crisis.electrodo_init(day)/s_rate;
        electrodo_end=patients{p}.crisis.electrodo_end(day)/s_rate;
        for i=1:length(header.channels)
            signal=data(i,:);
            window_size = s_rate*Time;
            overlap = 0;
            [Y,F,T,PSD] = spectrogram(signal,window_size,overlap,Frequencies,s_rate,'power','yaxis');
            delta=sum(PSD(1*10:3.5*10,:));
            theta=sum(PSD(3.5*10:7.4*10,:));
            alpha=sum(PSD(7.4*10:12.4*10,:));
            beta=sum(PSD(12.4*10:24*10,:));
            gamma=sum(PSD(24*10:97*10,:));
            matrix=[delta;theta;alpha;beta;gamma];
            %matrix2=matrix-repmat(mean(matrix,2),1,b);
            tiempo=0:Time:length(matrix)*Time-Time;
            count1=1;
            count2=1;
            for j=1:length(tiempo)
               if tiempo(j)<electrodo_init || tiempo(j)>electrodo_end
                  data_no_crisis(:,count1)=matrix(:,j);
                  count1=count1+1;
               else
                  data_crisis(:,count2)=matrix(:,j);
                  count2=count2+1;
               end
            end
            [aux1 aux2]=size(data_no_crisis);
            data_no_crisis2=data_no_crisis-repmat(mean(data_no_crisis,2),1,aux2);
            [V,D]=eig(data_no_crisis2*data_no_crisis2');%calculos los autovectores sin crisis
            [aux1 aux2]=size(data_crisis);
            data2=(D^-0.5)*V'*(data_crisis-repmat(mean(data_no_crisis,2),1,aux2));%proyecto los vectores con crisis en el espacio sin crisis
            [V_crisis,D_crisis]=eig(data2*data2'); % calculo autovectores con crisis
            eigenvector{p}(number,i,:,:)=V_crisis;%guardo los autovectores con crisis
            eig_crisis{p}(number,i,:)=D_crisis(find(D_crisis));%guardo los autovalores con crisis
            eigenvalues{p}(number,i,:,:)=RotationTest2(data2,100);%guardo el analisis de significancia
            max_eigenvalues=max(reshape(eigenvalues{p}(number,i,:,:),[5 100])');%calculo los valores máximos en las rotaciones de significancia
            sign=find(max_eigenvalues<reshape(eig_crisis{p}(number,i,:),[5 1])');%busco los autovalores significativos
            significance_position=zeros(5,1);
            significance_position(sign)=1;
            significance{p}(number,i,:)=significance_position;
            sign_vector=V*(D^-0.5)*V_crisis;%+repmat(mean(data_no_crisis,2),1,5); %reproyecto los vectores en crisis al espacio sin crisis
            vector_projection{p}(number,i,:,:)=sign_vector;            
        end
    end
end

nombre=strcat(directory,'eig_significance2')
save(nombre,'eigenvector','eig_crisis','eigenvalues','significance','vector_projection','-mat')