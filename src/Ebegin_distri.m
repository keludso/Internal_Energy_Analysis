function E = Ebegin_distri( Ebegin, Start_distr, Temp, Step )


ind = round(Temp/100);

if ind > 15 
    ind = 15;
end

 r_no = rand(1);
 k = 1;
 E_dist = Start_distr(ind,:);
 begin_st = Start_distr(k);
 k = k+1;
 while(begin_st<=r_no)
    
     begin_st = begin_st +edist(d);
     k = k+1;
        
 end
 
 E = k*Step;


end

