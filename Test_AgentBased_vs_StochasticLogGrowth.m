% This code is intended to check tumor growth and behaviour with 2
% different modelling frameworks: agent-based modelling, and stochastic
% logistic equation. Both of them should offer the same results; the point
% is checking which one works the fastest

clear all; close all; clc;


%% AGENT-BASED MODEL
%--------------------------------------------------------------------------

clear all; close all; clc;

tic 

% SET SPATIAL GRID
%--------------------------------------------------------------------------
N = 21;                 % Number of voxels per side. For a 3D system, total number of voxels will be N^3
x = linspace(1,N,N);    % Create arrays for each dimension (just to iterate easily along them)
y = x; 
z = 1; %z = y;
grid = meshgrid(x,y,z); % Set a spatial grid with specified dimensions

% SET NUMBER OF ALTERATIONS
%--------------------------------------------------------------------------
Alt = 6;

% SET TIME STEPS
%--------------------------------------------------------------------------
% We consider a year of time, evaluated each 10 hours
deltat = 10; % hours
T = 2*365*24;  % hours
Nstep = round(T/deltat);


% SET REQUIRED CELL LISTS
%--------------------------------------------------------------------------
Vs = cell(size(grid));  % Initialize vector state (with grid size, at each voxel contains 4 parameters) 
G = cell(size(grid));   % Initialize genotypes (each voxel will have a vector with genotype populations;
                        % genotypes will be defined by binary version of vector indexes)


% FILL ARRAYS AND CELL LISTS WITH INITIAL DATA
%--------------------------------------------------------------------------
R = 1;      % Resources
Nec = 0;    % Necrotic cells
P = 1e2;    % Cell population
K = 1e6;    % Limit population
for i = x
    for j = y
        for k = z
            Vs{i,j,k} = [R,Nec,0,K];    % Every position in vector state will contain 4 elements
            G{i,j,k} = zeros(1,2^Alt);  % Every position in genotype list will contain a vector of 2^Alterations elements
            %G{i,j,k}(1) = P;            % Initially, only healthy tissue will exist in a given voxel
        end
    end
end

G{11,11,1}(1) = P; % There will only be cells in central position
Vs{11,11,1}(3) = P;
% This replicates are made for coherent updating
Gnext = G;              
Vsnext = Vs;

% ASSOCIATED VARIABLES
%--------------------------------------------------------------------------
% Used to monitor system evolution
pop = zeros(N,N);
Neval = round(Nstep/1);
totpop = zeros(1,Neval);

% TIME EVOLUTION
%--------------------------------------------------------------------------
evalstep = 1;
for t = [1:Nstep]
    
    % Reinitialize variable to store population at certain time steps, as
    % it is going to be calculated before
    if mod(t,Nstep/Neval) == 0
        totpop(evalstep) = 0;
    end
    
    for i = x
        for j = y
            for k = z
                % 'e' stands for each genotype existing in a given voxel
                for e = [1:size(G{i,j,k},2)] 
                    if G{i,j,k}(e) > 0 % Only if there is population for the given genotype
                        % Reproduction event
                        Gnext{i,j,k}(e) = Gnext{i,j,k}(e) + rep_ab(G{i,j,k}(e),Vs{i,j,k}(3),Vs{i,j,k}(4));
                        % Death event
                        Gnext{i,j,k}(e) = Gnext{i,j,k}(e) - kill_ab(G{i,j,k}(e),Vs{i,j,k}(3),Vs{i,j,k}(4));
%                         % Migration event
%                         Pmig = probmig(migfit([i,j,k],Gnext,Vsnext,e));
%                         Pmig = Pmig*deltat;
%                         Gnext = mig(Gnext,Vsnext,[i,j,k],e,Pmig);

                    end
                end
                
                % Update vector state
                Vsnext{i,j,k}(3) = 0;
                for e = [1:size(Gnext{i,j,k},2)]
                    if G{i,j,k}(e) > 0 % Only if there is population for the given genotype
                        Vsnext{i,j,k}(3) = Vsnext{i,j,k}(3)+ Gnext{i,j,k}(e);
                    end
                end
                if mod(t,Nstep/Neval) == 0
                    totpop(evalstep) = totpop(evalstep)+Vsnext{i,j,k}(3);
                end
            end
        end
    end 
    
    % Update system
    G = Gnext;
    Vs = Vsnext;
    
    % Display some information to check how everything works
    if mod(t,Nstep/Neval) == 0
        evalstep = evalstep+1;
        disp(['Iteration nº ' num2str(t)])
        disp(['Total population of random voxel: ' num2str(Vsnext{11,11,1}(3))])
    end
end

toc

% Plots
%--------------------------------------------------------------------------

figure()
plot(([1:evalstep-1]*2)*10/24,totpop,'LineWidth',1)
xlabel('Time','FontSize',18)
ylabel('N cells','FontSize',18)
title(['Pop. growth with %_{prolif} = ' num2str(0.01) ', P_{kill} = ' num2str(0)],'FontSize',16)



%% STOCHASTIC LOGISTIC MODEL
%--------------------------------------------------------------------------

clear all; close all; clc;

tic 

% SET SPATIAL GRID
%--------------------------------------------------------------------------
N = 21;                 % Number of voxels per side. For a 3D system, total number of voxels will be N^3
x = linspace(1,N,N);    % Create arrays for each dimension (just to iterate easily along them)
y = x; 
z = 1; %z = y;
grid = meshgrid(x,y,z); % Set a spatial grid with specified dimensions

% SET NUMBER OF ALTERATIONS
%--------------------------------------------------------------------------
Alt = 6;

% SET TIME STEPS
%--------------------------------------------------------------------------
% We consider a year of time, evaluated each 10 hours
deltat = 10; % hours
T = 2*365*24;  % hours
Nstep = round(T/deltat);


% SET REQUIRED CELL LISTS
%--------------------------------------------------------------------------
Vs = cell(size(grid));  % Initialize vector state (with grid size, at each voxel contains 4 parameters) 
G = cell(size(grid));   % Initialize genotypes (each voxel will have a vector with genotype populations;
                        % genotypes will be defined by binary version of vector indexes)


% FILL ARRAYS AND CELL LISTS WITH INITIAL DATA
%--------------------------------------------------------------------------
R = 1;      % Resources
Nec = 0;    % Necrotic cells
P = 1e2;    % Cell population
K = 1e6;    % Limit population
for i = x
    for j = y
        for k = z
            Vs{i,j,k} = [R,Nec,0,K];    % Every position in vector state will contain 4 elements
            G{i,j,k} = zeros(1,2^Alt);  % Every position in genotype list will contain a vector of 2^Alterations elements
            %G{i,j,k}(1) = P;            % Initially, only healthy tissue will exist in a given voxel
        end
    end
end

G{11,11,1}(1) = P; % There will only be cells in central position
Vs{11,11,1}(3) = P;
% This replicates are made for coherent updating
Gnext = G;              
Vsnext = Vs;

% ASSOCIATED VARIABLES
%--------------------------------------------------------------------------
% Used to monitor system evolution
pop = zeros(N,N);
Neval = round(Nstep/2);
totpop = zeros(1,Neval);

% TIME EVOLUTION
%--------------------------------------------------------------------------
evalstep = 1;
for t = [1:Nstep]
    
    % Reinitialize variable to store population at certain time steps, as
    % it is going to be calculated before
    if mod(t,Nstep/Neval) == 0
        totpop(evalstep) = 0;
    end
    
    for i = x
        for j = y
            for k = z
                % 'e' stands for each genotype existing in a given voxel
                for e = [1:size(G{i,j,k},2)] 
                    if G{i,j,k}(e) > 0 % Only if there is population for the given genotype
                        % Reproduction event
                        Gnext{i,j,k}(e) = Gnext{i,j,k}(e) + rep_sl(G{i,j,k}(e),Vs{i,j,k}(3),Vs{i,j,k}(4));
%                         % Migration event
%                         Pmig = probmig(migfit([i,j,k],Gnext,Vsnext,e));
%                         Pmig = Pmig*deltat;
%                         Gnext = mig(Gnext,Vsnext,[i,j,k],e,Pmig);

                    end
                end
                
                % Update vector state
                Vsnext{i,j,k}(3) = 0;
                for e = [1:size(Gnext{i,j,k},2)]
                    if G{i,j,k}(e) > 0 % Only if there is population for the given genotype
                        Vsnext{i,j,k}(3) = Vsnext{i,j,k}(3)+ Gnext{i,j,k}(e);
                    end
                end
                if mod(t,Nstep/Neval) == 0
                    totpop(evalstep) = totpop(evalstep)+Vsnext{i,j,k}(3);
                end
            end
        end
    end 
    
    % Update system
    G = Gnext;
    Vs = Vsnext;
    
    % Display some information to check how everything works
    if mod(t,Nstep/Neval) == 0
        evalstep = evalstep+1;
        disp(['Iteration nº ' num2str(t)])
        disp(['Total population of random voxel: ' num2str(Vsnext{11,11,1}(3))])
    end
end

toc

% Plots
%--------------------------------------------------------------------------

figure()
plot(([1:evalstep-1]*2)*10/24,totpop,'LineWidth',1)
xlabel('Time','FontSize',18)
ylabel('N cells','FontSize',18)
title(['Pop. growth with t_{dup} = ' num2str(48) ' h'],'FontSize',16)



%% FUNCTIONS
%--------------------------------------------------------------------------

%% REPRODUCTION

% Everytime a reproduction event takes place, the chosen genotype must 
% increase its population, depending on a probability calculated from
% current population and limit size

% Reproduction of a given genotype for agent-based model
%--------------------------------------------------------------------------
function Popnew = rep_ab(Popgen,Pop,K)
    
    Prep = (1-Pop/K);
    Popnew = 0;
    
    % Not all cells are proliferative
    Percent_prolif = 0.02;
    
    if Prep > 0
        
%         % 1st way: using vectors
%         Rands = rand(1,round(Percent_prolif*Popgen));
%         Grow = find(Rands < Prep);
%         Popnew = Popnew + length(Grow);
        
        % 2nd way: using loops
        for i = [1:round(Percent_prolif*Popgen)]
            r = rand;
            if r < Prep
                Popnew = Popnew+1;
            end
        end

    end
    
end

% Reproduction of a given genotype for stochastic-logistic model
%--------------------------------------------------------------------------
function Popnew = rep_sl(Popgen,Pop,K)
    
    deltat = 10/24; r = 1/36;

    Ct = @(C,deltat,K,r) ((K*C*exp(r*deltat))/(K+C*(exp(r*deltat)-1)));
    Popnew = Ct(Popgen,deltat,K,r);
    Popnew = normrnd(Popnew,2e3*Popnew/K) - Popgen; % This way, the bigger variability is achieved in the end
    %Popnew = normrnd(Popnew,2e-3*K/Popnew) - Popgen; % This way, the bigger variability is achieved at the start
    
 
    Popnew = round(Popnew);

    
end


%% DEATH

% Everytime a death event takes place, the chosen genotype must decrease
% its population, depending on a uniform probability which is set constant

% Death of a given genotype for agent-based model
%--------------------------------------------------------------------------
function Popkill = kill_ab(Popgen,Pop,K)
    
    Pkill = 0.01;
    Popkill = 0;
    
    if Pkill > 0
        
        % 1st way: using vectors
        Rands = rand(1,Popgen);
        Kill = find(Rands < Pkill);
        Popkill = Popkill + length(Kill);
        
%         % 2nd way: using loops
%         for i = [1:Popgen]
%             r = rand;
%             if r < Pkill
%                 Popkill = Popkill+1;
%             end
%         end

    end
    
end


% Death of a given genotype for stochastic-logistic model (NOT USED)
%--------------------------------------------------------------------------
function Popkill = kill_sl(Popgen,Pop,K)
    
    Pkill = 0.01;
    Popkill = 0;
    
    if Pkill > 0
        
        % 1st way: using vectors
        Rands = rand(1,Popgen);
        Kill = find(Rands < Pkill);
        Popkill = Popkill + length(Kill);
        
%         % 2nd way: using loops
%         for i = [1:Popgen]
%             r = rand;
%             if r < Pkill
%                 Popkill = Popkill+1;
%             end
%         end

    end
    
end


