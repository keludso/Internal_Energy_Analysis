classdef MonteCarlo
    % A library for constants and functions 
    
    properties
        
        %% Input Parameters 
        Step =0;            % Grain Size
        TLIM = 0;
        Emax =0;
        Nmax =0;
        ALNDEN = zeros(1,4000);
        h = 1.9863e-23;
        Temp, TempI = 0;
        Ebegin = 0;
        Nstart = 0;
        DC = zeros(2,8);
        ITYPE   = [1,1];
        MW_power = 0;
        SIGMA0 = 0;
        Time_res =0;
        MW_on = 0;
        MW_Freq = 0;
        MW_absorbed = 0;
        Traj_react = zeros(2,4000);
        Edistribution= 0;
        
    end
    
end

