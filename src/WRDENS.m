function DENS = WRDENS(E, ALNDEN)

%     INTERPOLATE BETWEEN VALUES OF ALNDEN TO OBTAIN STATE DENSITY
%     BY USING QUADRATIC FITTING PROCEDURE
%     CALCULATE NEAREST ELEMENT IN ARRAY AND USE ADJACENT ELEMENTS

if E <= 0
   DENS = exp(ALNDEN(1));
   % WArig()
   return
end

IE = 1 + round(E/25);			% Nearest integer
   
if IE < 2 
   IE = 2;
end
   
if IE >= length(ALNDEN)
   IE = length(ALNDEN) -1;
end
       
X = E/25 - (IE - 1);

Y1 = exp(ALNDEN(IE-1));
Y2 = exp(ALNDEN(IE));
Y3 = exp(ALNDEN(IE+1));

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


if (Y <= 0)
    Y = 0;
    
end

DENS = Y;

end

