function [TAU, P] = Stochastic_calc(RKC1,RKC2, E,EC, ALNRKE, DC, Collup, Cnorm,...
                            ALNDEN, ICOLL, Temp, MW_on, MW_Freq, MW_power, Step)

  
%% Calculating reaction rates, timestep and probability of each reaction
%
% FIND TAU FOR CASE WHEN NO OPTICAL PUMPING (IFLAG1 = 0)
%

RKE = RRKM_rate(E, EC, ALNRKE, ALNDEN);

TARGET = log(1/rand(1));   % for gillespies

SUM = RKC1 + RKC2 + RKE;    % Total rate (no pump)
TAU = TARGET/SUM;         
P(1) = 0.0;     % Optical Pump
P(2) = 0.0;     % Optical Pump
      
if MW_on == 1    % Microwave 
     
      RXY = rand(1);
      Rnterp = RNTERP(E,ICOLL, Step,Collup);
     if (RXY < Rnterp || E <= 0)    	% An up-step was chosen
        IFLAG2 = 1;                % Collisional Up-step 
     else
        IFLAG2 =0;
     end
     STERAR = E_trans(E, DC, ICOLL, Step, Temp, ALNDEN, IFLAG2, Cnorm, Collup);
     rho = WRDENS(E + STERAR, ALNDEN);
     rho2 = WRDENS(E, ALNDEN);
     Ratio = rho/rho2;
     PUMPD = MW_power*0.318*(0.12/((STERAR - MW_Freq)^2 + 0.12^2));
     PUMPU = Ratio*PUMPD;
     SUM1 = SUM + PUMPU + PUMPD;
     TAU = TARGET/SUM1;
     P(1) = PUMPU*TAU;
     P(2) = PUMPD*TAU;
       
end

     P(3) = RKC1*TAU;        % Collision type 1
     P(4) = RKE*TAU ;       % Unimolecular rxn
     P(5) = RKC2*TAU;        % Collision type 2

     SUM = sum(P);
     P(1) = P(1)/SUM;

     for i = 2:5
       P(i) = P(i-1) + P(i)/SUM;    % P's compared later with random
     end                           %   numbers: to choose channel


end

