function CEI=CumulativeEpileptogenicityIndex(EI,v_bias)

    meanEI=mean(EI);
    
    for i=1:length(EI)
        CEI(i)=sum(EI(1:i))-i*meanEI-i*v_bias;
    end

end