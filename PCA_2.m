clear all
close all
directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/';
cd(directory)
data_base

directory='/home/melisa/Escritorio/Melisa/Proyecto_Epilepsia/PCA2/';
Frequencies=0:0.1:100;
bands=[1,3.5,7.4,12.4,24,97];
Time=1;
crisis={[2,3,5,6,7],[2,3],[1,3,4]};

for p=1:3
    last_end=1;
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
            data_no_crisis2=(data_no_crisis-repmat(mean(data_no_crisis,2),1,aux2));
            [V,D]=eig((data_no_crisis2*data_no_crisis2')/(length(data_no_crisis)-1));
            [aux1 aux2]=size(data_crisis);
            data2=(D^-0.5)*V'*(data_crisis-repmat(mean(data_no_crisis,2),1,aux2));
            data3=(D^-0.5)*V'*(data_no_crisis-repmat(mean(data_no_crisis,2),1,size(data_no_crisis,2)));
            [V_crisis,D_crisis]=eig((data2*data2')/(length(data_crisis)-1));
            %[eigenvalues_test eigenvectors_test]=RotationTest2(data2,1000);
            eigenvector{p}(number,i,:,:)=V_crisis;
            eig_crisis{p}(number,i,:)=D_crisis(find(D_crisis));
            eigenvalues{p}(number,i,:,:)=RotationTest2(data2,100);
            eigenvalues2{p}(number,i,:,:)=RotationTest3(data3,length(data2),100);          
            
            max_eigenvalues=max(reshape(eigenvalues{p}(number,i,:,:),[5 100])');%calculo los valores máximos en las rotaciones de significancia
            sign=find(max_eigenvalues<reshape(eig_crisis{p}(number,i,:),[5 1])');%busco los autovalores significativos
            
            max_eigenvalues2=max(reshape(eigenvalues2{p}(number,i,:,:),[5 100])');
            sign2=find(max_eigenvalues2<reshape(eig_crisis{p}(number,i,:),[5 1])');
            
            significance_position=zeros(5,1);
            significance_position(sign)=1;
            
            significance_position2=zeros(5,1);
            significance_position2(sign2)=1;
            
            significance{p}(number,i,:)=significance_position;
            significance2{p}(number,i,:)=significance_position2;

            sign_vector=V*(D^-0.5)*V_crisis;%+repmat(mean(data_no_crisis,2),1,5); %reproyecto los vectores en crisis al espacio sin crisis
            vector_projection{p}(number,i,:,:)=sign_vector;              
            figure(1)
            plot(reshape(eigenvalues{p}(number,i,:,:),[5 100]))
            hold on
            plot(reshape(eigenvalues2{p}(number,i,:,:),[5 100]))
            plot(reshape(eig_crisis{p}(number,i,:),[5 1]),'Linewidth',6)
            xlabel('Rank')
            ylabel('Eigenvalue')
            grid on
            box on
            title(strcat(patients{p}.name,'.Crisis=',int2str(crisis{p}(number)),'.Electrode=',header.channels(i,5:8)))
            name=strcat(directory,patients{p}.name,'.Crisis=',int2str(crisis{p}(number)),'.Electrode=',header.channels(i,5:8));
            saveas(1,name,'png')
            close (1)
        end
    end
end
nombre=strcat(directory,'eig_significance2')
save(nombre,'eigenvector','eig_crisis','eigenvalues','eigenvalues2','significance','significance2','vector_projection','-mat')

