%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo Simulation of Master equation to include microwave effects
% on the rate of the reaction.
%
% Author - Kelvin Dsouza, Dr Daryoosh Vashaee 
%
% Comment - The code is is converted from the MULTIWELL code in fortran
%  which results in calculating the unimolecular dissociation by Dr Barker
%
% The code begins with setting the parameter s for the reaction to occur. 
% 1. Number of trials - Number of trajectories in the path of the reaction 
% 2. Max Energy - Energy limit for simulation default- 100000
% 3. Grain Size - Energy grain size that would be used in the simulation. 
%               smaller the grain size better resolution 
% 4. Time Limit - The max time for simulation 
% 5. Microwave ON/OFF - Choosing the microwave settings 
% 6. Microwave Power - The amount of Fluence deposited on the sample during
%                       simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; tic

%% Input Configurations dialog box

prompt = {'File-name', 'Number of Trials','MAX Energy (cm-1)','Grid Size (cm-1):','Time Limit', ...
    'Microwave ON/OFF (1 or 0)', 'Microwave Power (J/cm2/s)', 'Starting Energy (cm-1)','Temperature (K)' ,'Frame size','Vibration data file (.txt)','Reaction Parameters (.txt)'};
title = 'Energy Distribution parameters';
dims = [1 60];
definput = {'name_me', '100','100000','25', '1e-6', '0','10', '40000','400', '100','vib.txt','reaction_parameters.txt' };
answer = inputdlg(prompt,title,dims,definput);

% Maximum energy and step size 

%% Global Variables 
% Assigning the Monte Carlo Class

cd src
Monte = MonteCarlo;

Monte.Nmax = str2double(answer(3));
Monte.Step = str2double(answer(4));
Monte.Emax = Monte.Nmax/Monte.Step;  % Maximum energy = Nmax*Step
Energy = 1:Monte.Step:Monte.Nmax;
Monte.Time_res = str2double(answer(10));

%% SET NUMBER OF TRIALS 
%
%  Nstart = NO. OF TRAJECTORIES  
%                         	     

Monte.Nstart = str2double(answer(2)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Calculating the density of states   %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[density, sum] = calculate_density(Monte.Emax,Monte.Step, answer{11});
Monte.ALNDEN = log(density);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET THE EBEGIN, FLAGS
% EBEGIN  - Starting energy point for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Monte.Ebegin   = str2double(answer(8));

%% SETTING TEMPERATURE 
%   TEMPIV = INTITAL VIBRATIONAL TEMPERATURE
%   TEMPI  = INITIAL TRANSLATIONAL TEMPERATURE
%     To simulate shock tube experiments, set TEMPIV equal to the temperature 
%     prior to the shock and set TEMPI equal to the initial (prior to relaxation) 
%     translational temperature.
%

Monte.TempI  = str2double(answer(9));
Monte.Temp   = Monte.TempI; % Temperature 


%% Energy Transfer Parameters 
%
%	Energy Transfer parameters
% DC Stores the coefficients to be used in the model

Monte.ITYPE   = [1,1];
Monte.DC(1,:) = [35.2, 0.0383, -1.18e-07, 1.5e-3,20000,0,0,0];
Monte.DC(2,:) = [28.4, 0.00521, -0.738e-07, 0,0,0,0,0];

%% Microwave Parameters

%   Microwave absorption parameters
%
%  IFLAG1 DESIGNATES OPTICAL PUMPING TYPE:
%
%     IFLAG1=
%          0:  NO OPTICAL PUMPING
%          1:  CONSTANT OPTICAL PUMPING,
%                FLUENC = POWER (J.CM-2.SEC-1)
%          2:  EXPONENTIAL PULSE, A1 IS DECAY TIME CONSTANT,
%               FLUENC IS TOTAL FLUENCE(J.CM-2)
%          3:  SQUARE PULSE OF FLUENC AND DURATION = 1/A1
%
%     MW_power = LASER PHOTON ENERGY (CM-1)
%

Monte.MW_on      = str2double(answer(6));   %MW ON/OFF
Monte.MW_power   = str2double(answer(7));
Monte.MW_Freq    = 8;

%%     TLIM IS THE TIME LIMIT OF THE CALCULATION
%

Monte.TLIM   = str2double(answer(5));
Monte.SIGMA0   = 1e-17;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Pre Processing 
%    Calculating Collision Rate 
%    Arranging THERMAL DISTR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulating the Monte Carlo Code 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Monte = Monte_Simulate_parfor(Monte,answer{12}, answer{1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
