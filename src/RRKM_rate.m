function [RKE] = RRKM_rate(E, EC, ALNRKE, ALNDEN)

RKE = 0;


if E > EC
    RKE = 0;
    EE = E - EC;

    if EE < 0 
        return
    else 
        IE = 1 + round(E/25);			% Nearest integer

        if IE < 2 
            IE = 2;
        elseif IE >= length(ALNDEN)
            IE = length(ALNDEN) -1;
        end

        X = E/25 - (IE - 1);
        Y1 = ALNRKE(1,IE-1);
        Y2 = ALNRKE(1,IE);
        Y3 = ALNRKE(1,IE+1);

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

        RKE = exp(Y);

    end
end

end

