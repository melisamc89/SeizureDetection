function  EI=EpileptogenicityIndex(F,step1,T,step2,PSD)

    count_t=0;
    t=0;
    while t < T(end)
        count_t=count_t+1;
        theta(count_t)=0;
        alpha(count_t)=0;
        beta(count_t)=0;
        gamma(count_t)=0;
        count=0;
        f=0;
        while f < F(end)
            count=count+1;
            if f > 3.5 & f < 7.4
                theta(count_t)=theta(count_t)+PSD(count,count_t);
            end
            if f > 7.4 & f <12.4
                alpha(count_t)=alpha(count_t)+PSD(count,count_t);
            end
            if f > 12.4 & f < 24
                beta(count_t)=beta(count_t)+PSD(count,count_t);
            end
            if f > 24 & f < 97
                gamma(count_t)=gamma(count_t)+PSD(count,count_t);
            end
            f=f+step1;
        end
        EI(count_t)=(beta(count_t)+gamma(count_t))/(theta(count_t)+alpha(count_t));
        t=t+step2;
    end

   
end
