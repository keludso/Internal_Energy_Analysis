function [ALNRKE, RKCI1, RKCI2, EC] = calculate_rate(Monte, f, file_name)



%% SETTING COLLISION PARAMETERS 
% PAPRES:
%   Two collider types are always assumed to be present: "Parent" and "backgrond
%   gas", or "collider".  To simulate an experiment containing 1.0 torr
%   benzene diluted in 50 torr of neon, PAPRES = 1.0 and BGP = 50.  Even if
%   only 1% of the benzene present is to be excited initially, the program will
%   assume that the complete 1.0 torr is excited.  This is often no problem,
%   except that the calculated temperature rise will be much too large.
%   The way to "trick" the program into giving the right result is to
%   increase all heat capacities and heat capacity derivatives by a factor
%   equal to the inverse of the fraction of the parent gas that is excited.
%
%  PARENT-GAS PRESSURE(TORR), Cp HEAT CAPACITY (CAL/DEG.MOLE) AT 300 K,
%         AND DERIVATIVE OF HEAT CAPACITY OF COLLIDER AND OF PARENT
%

% Details :
% The file should be in a txt file format and the following structure 
% 
% Line 1 :    "<Title>"
% Line 2 : 'Empty' Description for the user  
% Line 3 : <Number of vibrational mode>  <Number of Rotational Modes>
% Line 4 : <Mode = 'vib'/'rot'  <Frequency/ Ie>
% Line 5 : <Mode = 'vib'/'rot'  <Frequency/ Ie>
%

%%%%%%%% Constant Defination %%%%%%%%%%%%%%%%%%%

cd ..
cd input_files

rk_file = fopen(file_name);
i = 1;
tline = fgetl(rk_file);
title = tline;
comment = fgetl(rk_file);
tline = fgetl(rk_file);
i = 1;
while i<=16
    
   tline = fgetl(rk_file); 
   Str_file= tline;
   splt = strsplit(Str_file,' ');
   LJ_Param(i) = str2double(splt{2});
   i = i+1;
   
end 

fclose(rk_file);
cd ..
cd src

%% Assigning the values 
PAPRES = LJ_Param(1);
SIG1   = LJ_Param(2);
WELL1  = LJ_Param(3);
AMASS1 = LJ_Param(4);
CVPI   = LJ_Param(5);
CVPS   = LJ_Param(6);

% 0.1	3.655	178.9	84.	5.00	0.00		! Krypton

BGP     = LJ_Param(7);
SIG2    = LJ_Param(8);
WELL2   = LJ_Param(9);
AMASS2  = LJ_Param(10);
CVI     = LJ_Param(11);
CVS     = LJ_Param(12);

XSECT = 0.5*(SIG1 + SIG2);
WELL = sqrt(WELL1*WELL2);
RMASS = AMASS1*AMASS2/(AMASS1 + AMASS2);
BGP = BGP*9.66e18 /Monte.Temp;
PAPRES = PAPRES*9.66e18/Monte.Temp;
CV = (BGP*(CVI+(Monte.Temp - 300)*CVS) + PAPRES*(CVI + (Monte.Temp - 300)*CVPS))*5.808e-25;

waitbar(0.3,f,'Initialization.... Leonnard-Jones Potential');


% Reaction rate 1 
RKCI1 = (8.09e-10)*PAPRES*(sqrt(Monte.Temp/1000))*(sqrt(40/AMASS1))*((SIG1/5)^2)/(0.636 + (0.246*log(Monte.Temp/WELL1)));
% Reaction RAte 2 
RKCI2 = (8.09e-10)*BGP*(sqrt(Monte.Temp/1000))*(sqrt(20/RMASS))*((XSECT/5)^2)/(0.636+ 0.246*log(Monte.Temp/WELL));


%% RRKM parameters for UNIMOLECULAR DISSOCIATION

%     IRATE = 0 : INVERSE LAPLACE RRKM, 
%           = 1 : RRKM.  Input 300 element double-array of ln(sum of states) in ALNRKE
%     AVENU = LOG10(AFACTOR)
%     EC = critical energy (cm-1)
%     ROTDGN =  rotation factor*path degeneracy
%
% K(E) = A*density(E-Ecrit)/density(E)

IRATE   = LJ_Param(13);
AVENU   = LJ_Param(14);
EC      = LJ_Param(15);
ROTDGN  = LJ_Param(16);

ECRIT = Monte.Nmax - EC;

for j=1:length(Monte.ALNDEN)-2

       waitbar(j/length(Monte.ALNDEN),f,'Initialization.... RRKM rate');

       EE = (j-1)*Monte.Step;
       E  = EC + EE;

       DENS = WRDENS(E, Monte.ALNDEN); % Calculate density of states by interpolating the between values 
       DENOM = log(DENS);
       DENS = WRDENS(EE, Monte.ALNDEN);    % for Laplace Transform k(E)
       ANUM = log(DENS);
       ALNRKE(j) =  2.302585*AVENU + ANUM - DENOM;	% Inverse Laplace method      
end
        

delete(f);

end

