%% HETEROGENEITY IN GLIOBLASTOMA
%  Version 5.0
%
%  Use of indexes as binary coders for different genotypes has been
%  implemented -- as suggested by A. Martinez Rubio
%
%  Just check cell behaviour in 2D lattice, avoiding mutations and genotype
%  space interactions -- as suggested by Y. Azimzade
%
%  Use of binomial distribution to random sample newborn cells, instead of 
%  iterating along each individual -- as suggested by J. Jimenez

clear all; close all; clc;
tic


%% SETTING UP THE SYSTEM
%--------------------------------------------------------------------------

% SET SPATIAL GRID
N = 21;                 % Number of voxels per side. For a 3D system, total number of voxels will be N^3
x = linspace(1,N,N);    % Create arrays for each dimension (just to iterate easily along them)
y = x; 
z = 1; %z = y;
grid = meshgrid(x,y,z); % Set a spatial grid with specified dimensions

% SET NUMBER OF ALTERATIONS
Alt = 6;

% SET TIME STEPS
% We consider a year of time, evaluated each 10 hours
deltat = 10; % hours
T = 2*365*24;  % hours
%Nstep = T/deltat;
Nstep=500;


% SET REQUIRED CELL LISTS
Vs = cell(size(grid));  % Initialize vector state (with grid size, at each voxel contains 4 parameters) 
G = cell(size(grid));   % Initialize genotypes (each voxel will have a vector with genotype populations;
                        % genotypes will be defined by binary version of vector indexes)


% FILL ARRAYS AND CELL LISTS WITH INITIAL DATA
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

% ASSOCIATED VARIABLES
% Used to monitor system evolution
pop = zeros(N,N);
Neval = round(Nstep/2);
totpop = zeros(1,Neval);
totpop(1) = P;
totnec = zeros(1,Neval);
totnec(1) = 0;

% Duplicate cell list variables, to update system in a coherent way
Gnext = G;              
Vsnext = Vs;



%% LET THE SYSTEM EVOLVE
%--------------------------------------------------------------------------

evalstep = 1;
for t = [1:Nstep]
    
    for i = x
        for j = y
            for k = z
                % 'e' stands for each genotype existing in a given voxel
                for e = 1:size(G{i,j,k},2) 
                    if G{i,j,k}(e) > 0 % Only if there is population for the given genotype
                            
                        new = 0;
                        dead=0;

%                        Reproduction event
                        born = rep(G{i,j,k}(e),Vs{i,j,k}(3),Vs{i,j,k}(4),Vs{i,j,k}(2));
                        new = new + born;

%                         for l=1:Gnext{i,j,k}(e) % Iterate through each individual
%                        % Migration event
%                         [xmov,ymov] = mig(G,Vs,[i,j,k],e);
%                         new = new - 1;
%                         Gnext{xmov,ymov,k}(e) = Gnext{xmov,ymov,k}(e)+1;
%                         Vsnext{xmov,ymov,k}(3) = Vsnext{xmov,ymov,k}(3)+1;
%                         end

                        % Migration event 2/3
                        [Gnext,Vsnext,out]=mig2(G,Vs,[i,j,k],e,Gnext,Vsnext);
                        new=new-out;

                        % Death event
                        dead = death(G{i,j,k}(e),Vs{i,j,k}(3),Vs{i,j,k}(4),Vs{i,j,k}(2));
                        if dead < new
                            new = new - dead;
                        else
                            new = 0;
                            dead = new;
                        end


                        % Mutation event

                        % Updating
                        Vsnext{i,j,k}(3) = Vsnext{i,j,k}(3)+new;
                        Vsnext{i,j,k}(2) = Vsnext{i,j,k}(2)+dead;
                        Gnext{i,j,k}(e) = Gnext{i,j,k}(e)+new;
                      
                    end
                end
            end
        end
    end
    
    G = Gnext;
    Vs = Vsnext;
    
    if mod(t,Nstep/Neval) == 0
        totpop(evalstep+1) = 0;
        totnec(evalstep+1) = 0;
        for i = x
            for j = y
                for k = z
                    if Vsnext{i,j,k}(3) > 0
                        totpop(evalstep+1) = totpop(evalstep+1) + Vsnext{i,j,k}(3);
                        totnec(evalstep+1) = totnec(evalstep+1) + Vsnext{i,j,k}(2);
                    end
                end
            end
        end
        evalstep = evalstep+1;
         disp(['Iteration n∫ ' num2str(t)])
         disp(['Popgen: ' num2str(Vs{11,11,1}(3)) ', Totpop: ' num2str(totpop(evalstep)) ', Nec: ' num2str(Vs{11,11,1}(2))])
    end
end



%% PLOTS
%--------------------------------------------------------------------------

for i = x
    for j = y
        pop(i,j) = Vs{i,j,1}(3);
    end
end

toc

figure()
hold on
colormap gray
imagesc(pop)
axis tight
colorbar
xlabel('x')
ylabel('y')
title('Population per voxel')
hold off

figure()
teval = [0:1:length(totpop)-1];
hold on
plot(teval,totpop,'b-')
plot(teval,totpop+totnec,'r-')
xlabel('Time')
ylabel('Population')
title('Time evolution of total number of cells')
legend('P cells','P + N cells','Location','Best')
hold off




%% FUNCTIONS
%--------------------------------------------------------------------------

%% MUTATION

% Everytime a mutation event takes place, a new genotype will appear. The
% alteration of this new genotype will be chosen randomly. Only point
% alterations can happen; a jump of 2 or more alterations from a given
% genotype will never occur

% Mutation of a given genotype
%--------------------------------------------------------------------------




%% REPRODUCTION

% Everytime a reproduction event takes place, the chosen genotype must 
% increase its population in a random sampled quantity from a binomial 
% distribution, with a given probability of reproduction, and as much 
% number of trials as cells are attempting to reproduce

% Reproduction of a given genotype
%--------------------------------------------------------------------------
function born = rep(Popgen,Poptot,K,Nec)
    Prep = 0.01*(1-(Poptot+Nec)/K);
    born = binornd(Popgen,Prep);
end



%% MIGRATION

% Everytime a migration event takes place, the chosen genotype's population
% must lose an individual with a given probability. The individual will be 
% relocated on a nearby voxel. The voxel to be relocated will be chosen 
% randomly from the neighbourhood, depending on how crowded it is (the 
% emptier a voxel is, the higher the probability to be chosen)

% Migration of a cell for a given genotype
%--------------------------------------------------------------------------
function [xmov,ymov] = mig(gen,state,voxel,choice)

    x = voxel(1);
    y = voxel(2);
    z = voxel(3);
    
    Popgen = gen{x,y,z}(choice);
    Poptot = state{x,y,z}(3);
    K = state{x,y,z}(4);
    Nec = state{x,y,z}(2);
    
    Pmig = 1e-2;
    Lap = zeros(4,3);
    i = 0;
    for movi = [-1:1:1]
        for movj = [-1:1:1]
            xmov = x+movi;
            ymov = y+movj;
            
            if xmov < 21 && xmov > 0 && ymov < 21 && ymov > 0 && abs(movi)+abs(movj) ~= 0 && abs(movi)+abs(movj) ~= 2
                i = i + 1;
                lap = Pmig*(Poptot - state{xmov,ymov,z}(3));
                if lap < 0
                    lap = 0;
                end
                Lap(i,:) = [xmov,ymov,lap];
            end
        end
    end
    
    Lap = Lap(all(Lap,2),:);
    
    Lap(:,3) = Lap(:,3)/Poptot;
    
    probstay = 1-sum(Lap(:,3));
    
    Lap(length(Lap(:,3))+1,:) = [x,y,probstay];
    
    a = 1:length(Lap(:,1));
    
    select = a( sum( (rand(1) >= cumsum(Lap(:,3)./sum(Lap(:,3))))) + 1);
    
    xmov = Lap(select,1);
    ymov = Lap(select,2);

end



%% DEATH

% Everytime a death event takes place, the chosen genotype's population
% should be reduced in a random sampled quantity from a binomial 
% distribution, with a given probability of death, and as much 
% number of trials as cells are attempting to carry on apoptosis

% Death of a given genotype
%--------------------------------------------------------------------------
function dead = death(Popgen,Poptot,K,Nec)
    Pkill = 1e-4;%*((Poptot+Nec)/K);
    dead = binornd(Popgen,Pkill);
end


%% MIGRATION 2
%
% Select number of individuals according to pmig and then assign
% destination according to local laplacian
%
%--------------------------------------------------------------------------
function [gennext,statenext,migrants] = mig2(gen,state,voxel,choice,gennext,statenext)

    x = voxel(1);
    y = voxel(2);
    z = voxel(3);
    
    Popgen = gen{x,y,z}(choice);
    Poptot = state{x,y,z}(3);
    K = state{x,y,z}(4);
    Nec = state{x,y,z}(2);
    
    Pmig = 1e-2;
    Lap = zeros(4,3);
    i = 0;
    for movi = -1:1
        for movj = -1:1
            xmov = x+movi;
            ymov = y+movj;
            
            if xmov < 22 && xmov > 0 && ymov < 22 && ymov > 0 && abs(movi)+abs(movj) ~= 0 && abs(movi)+abs(movj) ~= 2
                i = i + 1;
                lap = (Poptot - state{xmov,ymov,z}(3));
                if lap < 0
                    lap = 0;
                end
                Lap(i,:) = [xmov,ymov,lap];
            end
        end
    end
    
    Lap(:,3) = Lap(:,3)/sum(Lap(:,3));
    
    a = 1:length(Lap(:,1));
    
    % Sample number of migrants and assign destination according to
    % probabilities
    
    migrants = binornd(Popgen,Pmig);
    
    if migrants > 0
    
    for i=1:migrants 
    
    select = a( sum( (rand(1) >= cumsum(Lap(:,3)./sum(Lap(:,3))))) + 1);
    
    xmov = Lap(select,1);
    ymov = Lap(select,2);
    
    gennext{xmov,ymov,z}(choice) = gennext{xmov,ymov,z}(choice)+1;
    statenext{xmov,ymov,z}(3) = statenext{xmov,ymov,z}(3)+1;
    
    
    end
    
    end

end


%% MIGRATION 3
%
% Compute local laplacian and sample individual according to multinomial
% distribution (with 4 possible outcomes - 4 neighbours)
% 
%--------------------------------------------------------------------------



function [gennext,statenext,out] = mig3(gen,state,voxel,choice,gennext,statenext)

    x = voxel(1);
    y = voxel(2);
    z = voxel(3);
    
    Popgen = gen{x,y,z}(choice);
    Poptot = state{x,y,z}(3);
    K = state{x,y,z}(4);
    Nec = state{x,y,z}(2);
    
    Pmig = 1e-2;
    Lap = zeros(4,3);
    i = 0;
    for movi = -1:1
        for movj = -1:1
            xmov = x+movi;
            ymov = y+movj;
            
            if xmov < 22 && xmov > 0 && ymov < 22 && ymov > 0 && abs(movi)+abs(movj) ~= 0 && abs(movi)+abs(movj) ~= 2
                i = i + 1;
                lap = Pmig*(Poptot - state{xmov,ymov,z}(3));
                if lap < 0
                    lap = 0;
                end
                Lap(i,:) = [xmov,ymov,lap];
            end
        end
    end
    
    Lap = Lap(all(Lap,2),:);
    
    Lap(:,3) = Lap(:,3)/Poptot;
      
    probstay = 1-sum(Lap(:,3));
    
    Lap(length(Lap(:,3))+1,:) = [x,y,probstay];
    
    % Sample number of migrants and assign destination according to
    % probabilities
    
    migrants = mnrnd(Popgen,Lap(:,3));
    
    Lap(:,4) = migrants;
    
    for i=1:length(Lap(:,3))-1
       
        xmov = Lap(i,1);
        ymov = Lap(i,2);
    
        gennext{xmov,ymov,z}(choice) = gennext{xmov,ymov,z}(choice)+Lap(i,4);
        statenext{xmov,ymov,z}(3) = statenext{xmov,ymov,z}(3)+Lap(i,4);
    
    end
    
    out=sum(migrants(1:length(Lap(:,3))-1));
end
