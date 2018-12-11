function RNTERP  = RNTERP(E,ICOLL, Step, XX)
%
%     Interpolate "double array" XX(2,500) by quadratic method
%

if E <= 0
   RNTERP = XX(ICOLL,1);
   return
end

IE = 1 + round(E/Step)	;		%Nearest integer

if IE < 2 
    IE = 2;
elseif IE > 3999
    IE = 3999;
end

X = E/Step - (IE - 1);

Y1 = XX(ICOLL,IE-1);
Y2 = XX(ICOLL,IE);
Y3 = XX(ICOLL,IE+1);

%
%     LINEAR/QUADRATIC INTERPOLATION, RETURNING VALUE OF Y AT VALUE OF X
%       Y1 AT X = -1
%       Y2 AT X =  0
%       Y3 AT X = +1
%
%	IF Y1=Y2, OR Y2=Y3, USE LINEAR TO PREVENT "OVERSHOOT";
%		 OTHERWISE, USE QUADRATIC 
%       QUADRATIC FORM:  Y = A*X**2 + B*X + C
%

if Y1~=Y2 && Y3 ~= Y2		% Quadratic interpolation
   C = Y2;
   B = 0.5*(Y3-Y1);
   A = 0.5*(Y1+Y3) - Y2;
   Y = (A*X + B)*X + C;  
elseif Y1 == Y2 && X <= 0	% Linear interpolation
   Y = Y2;
else
   Y = Y2 + X*(Y3-Y2);
end

RNTERP = Y;


end

