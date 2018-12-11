function [Dens, Ssum] = calculate_density(Emax,Step,file_name)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating the density and sum of states with the vibrational spectrum provided.
% 
%     S.E. STEIN AND B.S. RABINOVITCH, J.CHEM.PHYS. 58, 2438 (1973).
%   
%
%      Input: Emax = Number of energy grains 
%             Step = Step size of energy grain
%             File_name = Vibrational data file   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reading the vibrational Spectrum

% Details :
% The file should be in a txt file format and the following structure 
% 
% Line 1 :    "<Title>"
% Line 2 : 'Empty' Description for the user  
% Line 3 : <Number of vibrational mode>  <Number of Rotational Modes>
% Line 4 : <Mode = 'vib'/'rot'  <Frequency/ Ie>
% Line 5 : <Mode = 'vib'/'rot'  <Frequency/ Ie>
%
%
% ********** Read file **************************** 

cd ..
cd input_files

vib_lev = fopen(file_name);
i = 1;
tline = fgetl(vib_lev);
title = tline;
tline = fgetl(vib_lev);
tline = fgetl(vib_lev); 
Str_file= tline;
splt = strsplit(Str_file,' ');
N(i) = str2num(splt{1});
N_rots(i) = str2num(splt{2});

i = 1;
while i<=N
    
   tline = fgetl(vib_lev); 
   Str_file= tline;
   splt = strsplit(Str_file,' ');
   Mode(i) = str2num(splt{1});
   we(i) = str2num(splt{2});
   i = i+1;
   
end 

fclose(vib_lev);


cd ..
cd src
% ********** Close reading File **********
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **********  Density of States Calculation  **************************** 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Details : Density of states is calculated by the Stein Rabinovitch
% Algorithm. The Array T and AT are initialized to zero and the first trm
% is initialized to 1
% Refer to Robinson & Holbrook for equation

% Initializing Arrays 
T = zeros(1,Emax);
AT = zeros(1,Emax);
T(1) = 1;
AT(1) = 1;
Xe = 0;                  %Anharmonicity

% Running over each frequency 

for i = 1:N
    zpe = 0.5*we(i) + 0.25*Xe;                    % Zero point energy at v=0

    if Xe < 0
        w0 = we(i) + Xe;
        Rmax = -0.5*we(i)/Xe;    % MAximum bound level
    else
        Rmax = Emax*Step/we(i);
    end

    Nmax = round(Rmax);

    if Nmax > Emax
        Nmax = Emax;
    end

    for j = 1:Nmax
        vi = j +0.5;
        R_t = j*(we(i) + vi*Xe)*vi - zpe;   % ! state energy relative to zpe
        Ir(j) = 1+round(R_t/Step);
    end

    for j = 1: Nmax
        for k = Ir(j):Emax
            Karg  = k - Ir(j) + 1;
            if (Karg <= Emax)
                AT(k) = AT(k) + T(Karg);
            end
        end
    end

    for j = 1: Emax
        T(j) = AT(j) + T(j);
        AT(j) = 0;
    end

    % Calculating density
    % T(I) is the number of states in the Ith grain.

    Ssum(1) = T(1);                       %Energy at E = 0 (Energy at TOP of grain)
    Dens(1) = T(1)/Step;
    for j = 2:Emax
        Ssum(j)  = T(j) + Ssum(j-1);       % Sum of states (Energy at TOP of energy grain)
        Dens(j) = T(j)/Step;
    end
end
    
end

