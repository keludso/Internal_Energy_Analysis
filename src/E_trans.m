function STERAN = E_trans(E, DC, ICOLL, Step, TEMPI, ALNDEN, IFLAG2, Cnorm, Collup)

%     FROM RANDOM NUMBER RON AND ENERGY E COMPUTE A COLLISIONAL
%     STEP SIZE FROM THE NORMALIZED CDF CORRESPONDING TO
%     THE COLLISIONAL DOWN-STEP MODEL:
%
% 	IFLAG2=0 FOR DOWN-STEPS
% 	IFLAG2=1 FOR UP-STEPS
% 
% 	ICOLL=1  FOR PARENT-PARENT COLLISIONS
% 	     =2  FOR PARENT-COLLIDER COLLISIONS
% 
%        UPCOR = temperature correction for COLLUP (from TEMCOR)
%        CCOR  = temperature correction for CNORM  (from TEMCOR)
% 
% 	Step-size distribution given in FUNCTION PDOWN(E,EE,ICOLL,ITYPE,DC,TEMPI,SCALE)
% 
% 	Coefficients stored in array DC
% 	Probability of up-steps stored in COLLUP(2,500)
% 	Normalization factor (up + down) stored in CNORM(2,500)
%       Average energy transferred per collision stored in AEV(2,500)
% 
% 

Ediff = 1000;
error = 1.e-5;

deni = WRDENS(E, ALNDEN);			% DENI = Density at E
[SCALE, PROB] = PDOWN(E,E,DC, ICOLL, TEMPI);	% Evaluate to get initial SCALE
     
cnorme = RNTERP(E,ICOLL, Step,Cnorm);		%CNORME = Normalization at E
up= RNTERP(E,ICOLL,Step, Collup);

Tlast = PROB;				% "old" TERM
SUM = 0;
Test = 1;
EE = E;
r_no = rand(1);

if (IFLAG2 == 0) 			% Integrate over down steps
      
%******************** Down steps ***************

 ADOWN = r_no*cnorme*(1 - up);	%  Random number-selected integral of down-steps

    while SUM < ADOWN && EE  > 0 && Test > error
             
        H = (EE - Step * round(EE/Step));  	% Align to 25 cm-1 grain
        if H <= 0
            H = Step;
        elseif H > EE
            H = EE;
        end
        EE = EE - H;
        denj = WRDENS(EE, ALNDEN);
        if denj > 0.001
           [SCALE, PROB] = PDOWN(E,EE,DC, ICOLL, TEMPI);	% E > EE . % change E 
           TERM = PROB;
        else
           TERM = 0;
        end

        SPAN = 0.5 * H*(TERM + Tlast);

        if (SUM+SPAN) > ADOWN
            H = 1.01* (ADOWN - SUM)/(0.5* (TERM + Tlast));
            if H < Step
                H = Step;
            end 
            SPAN = 0.5 * H* (TERM + Tlast);  
        end
        SUM = SUM + SPAN;			% Down-step normalization integral
        Tlast = TERM;				% Old value for TERM

        if Tlast > 0.0 && abs(E-EE) > Ediff
           Test = abs(SPAN/SUM);
        end
        
    end

    if Tlast == 0
       EE = EE + H;
    end

        

else						% Up-collisions

% ******************** Up steps ***************

   ETOP = E + 10.*TEMPI/1.4388;			%  E + 10 kT
   AUP =  r_no*cnorme*up;			%  Random number-selected

   while SUM < AUP && EE<=ETOP && Test > error

      H = (EE - Step*round(EE/Step));  	% Align to 25 cm-1 grain
      if H <= 0
         H = Step;
      end
      
      EE = EE + H;
      denj = WRDENS(EE, ALNDEN);
          
      if denj >= 0.001
         [SCALE, PROB] = PDOWN(EE,E,DC, ICOLL, TEMPI);	% EE > E
         B = exp(-(EE-E)*1.4388/TEMPI);			% EE > E
         RATIO = (denj/deni);
         TERM = B*PROB*RATIO;
      else
         TERM = 0;
      end
          
      SPAN = 0.5*H*(TERM + Tlast);

      if (SUM+SPAN) > AUP 
         H = 1.01* (AUP - SUM)/(0.5*(TERM + Tlast));
         if ( H < Step) 
             H = Step;
         end
            SPAN = 0.5*H*(TERM + Tlast);
      end
          
      SUM = SUM + SPAN;			% Up-step normalization integral
      Tlast = TERM;			% Old value for TERM
      if (Tlast > 0 && abs(EE-E) > Ediff) 
          Test = abs(SPAN/SUM);
      end
          
     
  end
       
if (Tlast == 0) 
   EE = EE - H;
end
if EE < 0
   EE = 0;
end

end

STERAN = abs(E-EE);  
      
end
