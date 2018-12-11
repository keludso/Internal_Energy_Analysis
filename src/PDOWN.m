function [SCALE, PDOWN] = PDOWN(E, EE, DC, ICOLL, TEMPI)


%	PDOWN = step-size probability density distribution, not normalized,
%		for DOWN-steps:  E > EE
%	SCALE = characteristic "average" step size corresponding to PDOWN;
%		used in SUBROUTINE ANORML to determine integration step-sizes.
%
%	ICOLL = 1  for Parent + Parent
%	ICOLL = 2  for Parent + Collider



%  Biexponential Model                                                                                
%   Alpha1 = C(1) + [E*C(2) + E*E*C(3)]*T**C(8)                                                       
%   Alpha2 = C(5) + [E*C(6) + E*E*C(7)]*T**C(8)                                                       
%  Pdown = (1-C(4))*EXP(-(E-EE)/Alpha1) + C(4)*EXP(-(E-EE)/Alpha2)   
% 
% 


ALPHA1 = abs(DC(ICOLL,1) + (E*(DC(ICOLL,2) + E*DC(ICOLL,3)))*TEMPI^(DC(ICOLL,8)));
ALPHA2 = abs(DC(ICOLL,5) + (E*(DC(ICOLL,6) + E*DC(ICOLL,7)))*TEMPI^(DC(ICOLL,8)));

if ALPHA1 > 0.0
   TERM1 = exp(-abs(E-EE)/ALPHA1);
   TERM11 = TERM1/ALPHA1;
else
   TERM1 = 0;
   TERM11 = 0;
end


if ALPHA2 > 0
  TERM2 = exp(-abs(E-EE)/ALPHA2);
  TERM22 = TERM2/ALPHA2;
else
  TERM2 = 0;
  TERM22 = 0;
end

PDOWN = (1 - DC(ICOLL,4))*TERM1 + DC(ICOLL,4)*TERM2;

if PDOWN == 0
    SCALE = 0;
else   
    SCALE = PDOWN/((1 - DC(ICOLL,4))*TERM11 + DC(ICOLL,4)*TERM22);
end


if SCALE > max(ALPHA1,ALPHA2)
    SCALE = max(ALPHA1,ALPHA2);
elseif SCALE < min(ALPHA1,ALPHA2)
    SCALE = min(ALPHA1,ALPHA2);
end

      
end

