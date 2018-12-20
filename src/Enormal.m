function [Collup, Cnorm] = Enormal(ICOLL, Monte, Collup, Cnorm)


% c	Normalization for a down-step model given 
% c		in FUNCTION PDOWN(E,EE,ICOLL,ITYPE,DC,TEMPI,SCALE)
% c	
% c	where 	E = initial energy
% c		EE = final energy
% c		TEMPI is the temperature
% c		ICOLL designates parent (=1) or collider (=2)
% c		ITYPE = a parameter to define PDOWN
% c		DC are coefficients
% c		SCALE = characteristic (average) stepsize at energy E



error = 1.0D-4;
Ediff = 1000;	% Minimum integration interval (cm-1)
iter = 4;			% Iterations on up-steps
AED = zeros(1,4000);
AEV = zeros(1,4000);
Cn = zeros(1,4000);
%	Integrate over down steps
%		(temporarily store down-step integrals in COLLUP)

for i = 2:4000
      test = 1;
      AED(i) = 0;			% Storage of EE down integral
      Collup(ICOLL,i) =  0;	% Temporary storage of down-normalization
      E = (i-1)*Monte.Step;
      deni = WRDENS(E, Monte.ALNDEN);				% Density at E(I)
      
      if (deni > 0.001) 		%  ONLY IF STATE, OTHERWISE SKIP
          
        EE = E;
        [~, PROB] = PDOWN(E,E, Monte.DC, ICOLL, Monte.Temp);
        tlast = PROB;						% Initial SCALE and "old" TERM
        alast = 0;					% "old" AED

        while (test > error) && EE > 0	% ********** DOWN STEPS **********

             H = EE - 25*round(EE/25);  	% Align to 25 cm-1 grain
             if H <= 0
                 H = Monte.Step;
             elseif H > EE
                 H = EE;
             end
             
             EE = EE - H;
             denj = WRDENS(EE, Monte.ALNDEN);
          
            if denj >= 0.001
                [~, PROB] = PDOWN(E,EE, Monte.DC, ICOLL, Monte.Temp);	% E > EE
                term = PROB;
            else
                term = 0;
            end
            
            span = 0.5*H*(term + tlast);
            AED(i) = AED(i) + 0.5*H*((EE - E)*term + alast);	% Down-step delta-E integral
            Collup(ICOLL,i) = Collup(ICOLL,i) + span;		% Trapezoidal Rule
            if span > 0 && abs(E - EE)> Ediff
                test = abs(span/Collup(ICOLL,i));
            end
          
            alast = (EE - E)*term;					% "old" AED
            tlast = term;						% "old" TERM
          
        end
      end
      
end


AED(1) = 0;
Collup(ICOLL,1) =  0;
% 
% c	*************************************************************
% C	Now, Iterate over up-steps until normalization has converged
% 
% c	First approximation for normalization: down-integral + kT

for i = 1: 4000
    Cnorm(ICOLL,i) = Collup(ICOLL,i) + Monte.Temp/1.4388;
end

for k = 1:iter				% Iterations
    for i = 1:4000
       Cn(i) = Collup(ICOLL,i);		% Initialize by setting equal to down integral
       test = 1;
       cnorme = Cnorm(ICOLL,i);				% Normalization at E
       AEV(ICOLL,i) = AED(i);		% Initialize by setting equal to down integral

       E = (i-1)*Monte.Step;
       deni = WRDENS(E, Monte.ALNDEN);				% Density at E(I)

       if deni >=  0.001            %  PROCEED ONLY IF STATE PRESENT

           ETOP = E + 10*Monte.Temp/1.4388;				% E + 10 kT
           EE = E;
           [~, PROB] = PDOWN(E,E, Monte.DC, ICOLL, Monte.Temp);
           tlast = PROB	;					% Initial SCALE and "old" TERM
           alast = 0 ;   						% "old" AED

           while test > error && EE <= ETOP 	% ********** UP STEPS ***********
               H = (EE - Monte.Step*round(EE/Monte.Step)) ; 	% Align to 25 cm-1 grain
               
               if H <= 0
                  H = Monte.Step;
               end
               EE = EE + H ;
               denj = WRDENS(EE,Monte.ALNDEN);
          
               if (denj >= 0.001) 
                   [SCALE, PROB] = PDOWN(EE,E, Monte.DC, ICOLL, Monte.Temp);% EE > E
                   cnormee = RNTERP(EE,ICOLL, Monte.Step, Cnorm);		% Normalization at EE
                   B = exp(-(EE - E)*1.4388/Monte.Temp);			% EE > E
                   RATIO = cnorme*denj/(deni*cnormee);
                   term = B*PROB*RATIO;
               else
                   term = 0;
               end
          
               span = 0.5*H*(term + tlast);
               Cn(i) = Cn(i) + span;					% Trapezoidal Rule
               AEV(ICOLL,i) = AEV(ICOLL,i) + 0.5*H*((EE-E)*term + alast);
               if span > 0.0 && abs(E - EE) > Ediff
                   test = abs(span/Cn(i));
               end
          
               alast = (EE - E)*term;					% "old" value for AED
               tlast = term;						% "old" value of TERM
            end
       
       end
       
    end             % end energy steps

    for j = 1:4000
            Cnorm(ICOLL,j) = Cn(j);			%  Re-load temporary CN (normalization)

    end
end             % end iterations

% Now calculate final COLLUP array total delta-E: AEV

      for i = 1:4000
            if  Cnorm(ICOLL,i) > 0
                Collup(ICOLL,i) = 1 - Collup(ICOLL,i)/Cnorm(ICOLL,i);
                AEV(ICOLL,i) = AEV(ICOLL,i)/Cnorm(ICOLL,i);
            else
                Collup (ICOLL,i) = 0;
                AEV(ICOLL,i) = 0;
            end
            
      end
      

end

