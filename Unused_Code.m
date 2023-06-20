%Finding the error for limearly interpolaed values
ind1 = find(TP1.datenumb<t(j),1);
        ind2 = find(TP1.datenumb>t(j),1);
        %just do the average of two dates (initial and latest)
        err_upp(j) = sqrt(TP1.error_upp(ind1)^2 + TP1.error_upp(ind2)^2);
        error_upp(j) = sqrt(TP1.error_upp(ind1)^2+ err_upp(j)^2);
        err_low(j) = sqrt(TP1.error_down(ind1)^2 + TP1.error_down(ind2)^2);
        error_down(j) = sqrt(TP1.error_down(ind1)^2 + err_low(j)^2);


         %weighted error
        a = 0;b = 0;
        d1 = find(TP1.datenumb<t(j),1);
        dss = TP1.datenumb(d1);
        for i = 1:len
            dif = (TP1.datenumb(terms(i))-(t(j)-1+dif));
            avg_err_upp(i) = ((dif/days)*TP1.error_upp(terms(i)))^2;
            avg_err_low(i) = ((dif/days)*TP1.error_down(terms(i)))^2;
            a = a + avg_err_upp(i);
            b = b + avg_err_low(i);
            dss = dss+dif;
        end
        dif2 = t(j+1)-TP1.datenumb(terms(len));
        if dif2>0
            a = a+((dif2/days)*TP1.error_upp(terms(len)))^2;
            b = b+((dif2/days)*TP1.error_down(terms(len)))^2;
        end
        %calculate the difference between two dates
        err_upp(j) = sqrt(a)/dss;
        err_down(j) = sqrt(b)/dss;
     end