function Monte =  Monte_Simulate_parfor(Monte, r_file, comment)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% MonteCarlo Simulation 
% Simulation of MonteCarlo code 
% 
% Finding starting energy in the distribution %%
% Runs the gillespie algorithm to calculate the energy distribution based on
% reaction rate


%% Initialization 
 
f = waitbar(0,'Initialization..... ');
waitbar(0.3,f,'Calculating Density of states ..');
pause(0.1)

Step = Monte.Step;
TLIM = Monte.TLIM;
NET = 0;
MW_power  = Monte.MW_power;
traj = zeros(Monte.Nstart,1);
Ebegin = Monte.Ebegin;
ALNDEN = Monte.ALNDEN;
Emax = Monte.Emax;
Nmax = Monte.Nmax;
Ebegin = Monte.Ebegin;
Time_res = Monte.Time_res;
MW_on = Monte.MW_on;
MW_Freq = Monte.MW_Freq;

[ALNRKE, RKCI1, RKCI2, EC] = calculate_rate(Monte,f, r_file);
Collup = zeros(2,4000);
Cnorm = zeros(2,4000);
[Collup, Cnorm] = Enormal(1,Monte, Collup, Cnorm);
[Collup, Cnorm] = Enormal(2,Monte, Collup, Cnorm);
DC = Monte.DC;
NSTART =Monte.Nstart;
Temp = Monte.Temp;
NOTFIN = 0;
Coll_time =0;

f = waitbar(0,' Processor Allocation ....');
pause(.1)

NTRAJ_t =  cell(NSTART/100 ,1);
Coll_time_t = cell(NSTART/100 ,1);
Energy_DIST = cell(NSTART/100 ,1);
absorbed_avg =zeros(1,101);
emitted_avg =zeros(1,101);

%% Running the loop

for j = 1:NSTART/100 
    EDIST_par = cell(100 ,1);
    NTRAJ = zeros(1,100);
    Coll_time = zeros(1,100);  

   for i= 1:100
      ICOLL = 2; 
      tic;
      EDIST = zeros(Time_res+1,4000);
      Photon_absorbed = zeros(1,Time_res+1);
      Photon_emitted = zeros(1,Time_res+1);
      E = Ebegin;
      NET_A =0;
      NET_E =0;
      status = 0;
      T=0;

     
      while (status == 0)
        
          EE = E ;      % Energy for bookkeeping
          TT = T;       % Time for Book keeping
          if T >= TLIM % If time has excedded the limit next itteration 
            NOTFIN = NOTFIN + 1;
            break;
          end     
          [TAU, P] = Stochastic_calc(RKCI1,RKCI2, E,EC, ALNRKE, DC, Collup, Cnorm,...
                            ALNDEN, ICOLL, Temp, MW_on, MW_Freq, MW_power, Step);      
  
          if (TAU+ TT) > TLIM || isnan(TAU)
                TAU = 1.01*TLIM;
          end
          
     % Book keeping function
        
         T_start = round(Time_res*TT / TLIM) + 2; %    DIVIDE TIME INTO 100 INTERVALS FROM TIME = 0.0 TO TLIM
         T_start(T_start>Time_res+1) = Time_res+1;
         AARGU = (TT + TAU)/TLIM;
         T_finish = round(Time_res.*AARGU) + 2;
         T_finish(T_finish>Time_res+1) = Time_res+1;
         K = 1 + round(length(ALNDEN)*EE/Nmax);
     
        if K > length(ALNDEN)
            K = length(ALNDEN);
        end
        
        if TT == 0   %	AT T = 0 
            EDIST(1,K) = EDIST(1,K) + 1;
        end

        EDIST(T_start:T_finish,K) = EDIST(T_start:T_finish,K) + 1;
        Photon_absorbed(T_start:T_finish) = Photon_absorbed(T_start:T_finish) + NET_A;
        Photon_emitted(T_start:T_finish) = Photon_emitted(T_start:T_finish) - NET_E;
        
        if TAU == 0
            status = 0;
        elseif T > TLIM
            status = 1;
        else
            T = T+ TAU;
            [status, MW_absorbed, ICOLL, E, NET_A, NET_E] = monte_calc(P, DC, T, E, Collup, Cnorm, ...
                                MW_Freq, MW_on, Step, Temp, ALNDEN, MW_power);
        end
      end

    if T < TLIM 
        NTRAJ(i) = 1;
        Coll_time(i) = TT;       
    end
    
    EDIST_par{i}  = EDIST;
    Photon_absorbed_t{i} = Photon_absorbed;
    Photon_emitted_t{i} = Photon_emitted;
    MW_absorbed_t{i} = MW_absorbed;

  end

 avg_energy = 0;
 waitbar(100*j/NSTART,f,'Trials running... ');
 
for k = 1:100
    avg_energy = avg_energy + EDIST_par{k};
    absorbed_avg = absorbed_avg + Photon_absorbed_t{k};
    emitted_avg = emitted_avg + Photon_emitted_t{k};
    
end

Energy_DIST{j} = avg_energy;
NTRAJ_t{j} = NTRAJ;
Coll_time_t{j} = Coll_time;
    
end   
 
Traj = 1; 
avg_energy = zeros(101,4000);
Traj_react = zeros(1,101);

for k = 1:j
    avg_energy = avg_energy + Energy_DIST{k}; 
    NTRAJ_m = NTRAJ_t{k};
    Coll_time_m = Coll_time_t{k};
    for m = 1: 100      
        if NTRAJ_m(k) >= 1
            T_f  = round((100*Coll_time_m(k)/TLIM));
            Traj_react(T_f:end) = Traj;
            Traj = Traj + 1;
        end
    end  
end

Energy = 1:Step:100000; %max
Time = 0:TLIM/Time_res:TLIM;

delete(f);

% Storing the final output data
cd ..
cd data_file
file_name = strcat(string(comment),'.mat');
save (file_name);
cd ..



%% Saving to a file 

%data_write = fopen(file_name,'w');
%% Plotting the data 
figure;
for i = 1:101

    set(gcf,'color','w');
    subplot(2,2,3);
    title('Microwave Photons emitted');
    plot(Time(1:i), emitted_avg(1:i));
    xlabel('Time s');
    ylabel('Photons emitted');
    xlim([ 0 TLIM])
    pause(0.2);

    subplot(2,2,2);
    title('Trajectories Reacted');
    plot(Time(1:i), Traj_react(1:i));
    xlim([ 0 TLIM])
    xlabel('Time');
    ylabel('Trajectories reacted');

    subplot(2,2,1);
    title('Internal Energy Distribution');
    plot(Energy, avg_energy(i,:));
    xlabel('Energy cm^-^1');
    ylabel('Population');
    pause(0.2);

    subplot(2,2,4);
    title('Microwave Photons absorbed');
    plot(Time(1:i), absorbed_avg(1:i));
    xlabel('Time s');
    ylabel('Photons absorbed');
    xlim([ 0 TLIM])
    pause(0.2);
    
    


end

end

