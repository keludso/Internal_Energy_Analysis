function [status, MW_absorbed, ICOLL, E, NET_A, NET_B] = monte_calc(P, DC, T, E, Collup, Cnorm, ...
             MW_Freq, MW_on, Step, Temp, ALNDEN, MW_power)


%% This function calculates the reactiuon that will occur based on the
% probability of the reaction. 

ICOLL = 1;
NET_A = 0;
NET_B = 0;
MW_absorbed =0;
RX = rand(1);

if MW_on > 0
       if RX < P(1)                   % Photon absorption
            MW_jump = E_trans(E, DC, ICOLL, Step, Temp, ALNDEN, 1, Cnorm, Collup);
            E = E + MW_jump;
            NET_A = 1;
            status = 0; % Back to the start 
            return;
    
       elseif RX < P(2)        % Stimulated emission
           MW_jump = MW_Freq;
           E = E - MW_jump ;               % Photon stimulated emission
           if E < 0 
              E = 0;
           end
           NET_B = - 1;
           status =0;  % Back to the start 
           return;
       end
            
end

      if RX < P(3)           % Collision type 1
          RNTER = RNTERP(E,ICOLL, Step, Collup);
          RXY = rand(1);
          if (RXY < RNTER || E <= 0)    	% An up-step was chosen
             STERAR = E_trans(E, DC, ICOLL, Step, Temp, ALNDEN, 1, Cnorm, Collup);
             E = E + STERAR;
          else
            STERAR = E_trans(E, DC, ICOLL, Step, Temp, ALNDEN, 0, Cnorm, Collup);
            E = E - STERAR;
         end
         STERAR = E_trans(E, DC, ICOLL, Step, Temp, ALNDEN, 0, Cnorm, Collup);
         E = E - STERAR;
          
         status = 0;
         return;
          
      elseif RX < P(4)           % Unimolecular rxn
          
            MW_absorbed = MW_power * T;  
            status =1;
            return;
      
      elseif RX < P(5)           % Collision type 2
          ICOLL = 2;
          RNTER = RNTERP(E,ICOLL, Step,Collup);
          RXY = rand(1);          
         if (RXY < RNTER || E <= 0)    	% An up-step was chosen
               STERAR = E_trans(E, DC, ICOLL, Step, Temp, ALNDEN, 1, Cnorm, Collup);
               E = E + STERAR;
   
         else
          
            STERAR = E_trans(E, DC, ICOLL, Step, Temp, ALNDEN, 0, Cnorm, Collup);
            E = E - STERAR;
         end
         
         status = 0;
         return;           
     end
    
status =0;

end


