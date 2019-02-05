%% HETEROGENEITY IN GLIOBLASTOMA
%  Version 2.0

clear all; close all; clc;


%% SETTING UP THE SYSTEM
%--------------------------------------------------------------------------

% SET SPATIAL GRID
N = 5;                  % Number of voxels per side. For a 3D system, total number of voxels will be N^3
x = linspace(1,N,N);    % Create arrays for each dimension (just to iterate easily along them)
y = x; z = y;
grid = meshgrid(x,y,z); % Set a spatial grid with specified dimensions

% SET NUMBER OF ALTERATIONS
Alt = 6;

% SET REQUIRED CELL LISTS
Vs = cell(size(grid));  % Initialize vector state (with grid size, at each voxel contains 4 parameters) 
G = cell(size(grid));   % Initialize genotypes


% FILL ARRAYS AND CELL LISTS WITH INITIAL DATA
R = 1;    % Resources
N = 0;      % Necrotic cells
P = 1e2;    % Cell population
K = 1e4;    % Limit population
for i = x
    for j = y
        for k = z
            Vs{i,j,k} = [R,N,P,K];          % Every position in vector state will contain 4 elements
            G{i,j,k} = {zeros(1,Alt),P};    % Every position in genotype list will contain sets of 2 elements
                                            % 1st genotype will usually be healthy tissue
        end
    end
end




%% LET THE SYSTEM EVOLVE
%--------------------------------------------------------------------------

deltat = 10;
for t = [1:1000]
    
    % Duplicate cell list variables, to update system in a coherent way
    Gnext = G;              
    Vsnext = Vs;
    
    for i = x
        for j = y
            for k = z
                % 'e' stands for each genotype existing in a given voxel
                for e = [1:size(Gnext{i,j,k},1)] 
                    % Reproduction event
                    Prep = probfit(fitness([i,j,k],Gnext,Vsnext,e));
                    Prep = Prep*deltat;
                    Gnext = rep(Gnext,[i,j,k],e,Prep);
                    % Migration event
                    Pmig = probmig(migfit([i,j,k],Gnext,Vsnext,e));
                    Pmig = Pmig*deltat;
                    Gnext = mig(Gnext,[i,j,k],e,Pmig);
                    % Death event
                    Pkill = probkill([i,j,k],Gnext,Vsnext,e);
                    Pkill = Pkill*deltat;
                    Gnext = death(Gnext,[i,j,k],e,Pkill);
                    % Mutation event
                    Pmut = probmut([i,j,k],Gnext,Vsnext,e);
                    Pmut = Pmut*deltat;
                    Gnext = mut(Gnext,[i,j,k],e,Pmut);
                end
            end
        end
    end
    
    % Update population in state vector
    for i = x
        for j = y
            for k = z
                Vsnext{i,j,k}(3) = Gnext{i,j,k}{1,2};
            end
        end
    end   
    
    % Update system
    G = Gnext;
    Vs = Vsnext;
end



for i = x
    for j = y
        for k = z
            % Check if both values are the same (if not, IT'S A TRAP!!)
            %Gnext{i,j,k}{1,2}
            Vsnext{i,j,k}(3)
        end
    end
end


%% PLOTS
%--------------------------------------------------------------------------

figure()
hold on
for i = x
    for j = y
        %for k = z
            scatter3(i,j,N/2,1000,Vs{i,j,k}(3),'filled');
        %end
    end
end
colorbar
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
function gen = mut(gen,voxel,choice,Pmut)
    x = voxel(1);
    y = voxel(2);
    z = voxel(3);
    genotype = gen{x,y,z}{choice,1};

    mutated = 1;
    r = rand;
    if r < Pmut
        
        % Be sure of picking only a non-altered slot
        while genotype(mutated) == 1
            mutated = randi(length(genotype));
        end
        
        new_genotype = genotype;
        new_genotype(mutated) = 1;

        % Ancient genotype will lose an individual
        gen{x,y,z}{choice,2} = gen{x,y,z}{choice,2}-1;
        
        % New genotype will gain an individual
        %gen{x,y,z} = [gen{x,y,z};{new_genotype,1}];
        nG = size(gen{x,y,z},1);
        gen{x,y,z}{nG+1,1} = new_genotype;
        gen{x,y,z}{nG+1,2} = 1;
        disp('A MUTATION EVENT OCCURRED')
        x
        y
        z
    end
end


% Calculate probability of mutation
%--------------------------------------------------------------------------
function Pmut = probmut(voxel,gen,state,choice)
    x = voxel(1);
    y = voxel(2);
    z = voxel(3);
    
    Pmut = 1e-7;
    genotype = gen{x,y,z}{choice,1};

    % For each alteration the genotype already has, it is more likely for
    % it to mutate again
    for alt = [1,length(genotype)]
        if genotype(alt) == 1
            Pmut = Pmut + 1e-7;
        end
    end
    
    Pmut = Pmut;
end



%% REPRODUCTION

% Everytime a reproduction event takes place, the chosen genotype must 
% increase its population in 1, depending on a probability calculated from
% its fitness

% Reproduction of a given genotype
%--------------------------------------------------------------------------
function gen = rep(gen,voxel,choice,Prep)
    x = voxel(1);
    y = voxel(2);
    z = voxel(3);
    
    r = rand;
    if r < Prep
        gen{x,y,z}{choice,2} = gen{x,y,z}{choice,2}+1;
    end
end


% Calculate probability of reproduction
%--------------------------------------------------------------------------
function Prep = probfit(Hrep)
    Prep = 1/(1+exp(-Hrep+5));
    %Prep = 1/(1+exp(-Hrep));
end


% Calculate fitness
%--------------------------------------------------------------------------
function Hrep = fitness(voxel,gen,state,choice)
    % A 'voxel' contains all 3 coordinates needed to place it
    % 'gen' is the cell list containing all genotypes for each voxel
    % The 'state' is the state vector cell list for all voxels
    x = voxel(1);
    y = voxel(2);
    z = voxel(3);
    
    % Basic matrix for correlation between alterations (used to compute fitness) 
    M = [5.4,5,0,5,5,-5;5,4.6,0,0,0,-5;0,0,3.9,0,0,0;5,0,0,5.4,-5,5;5,0,0,-5,3.2,0;-5,-5,0,-5,0,4];
    M = (10/det(M)).*M;
    
    % This 
    nG = size(gen{x,y,z},1);
    genotype = gen{x,y,z}{choice,1};
    resources = state{x,y,z}(1)*state{x,y,z}(4);
    Hrep = genotype*M*genotype';
    for i = [1,nG]
        Gpop = gen{x,y,z}{i,2};
        Hrep = Hrep - 1/resources*Gpop;
    end

end


%% MIGRATION

% Everytime a migration event takes place, the chosen genotype's population
% must lose an individual with a given probability. The individual will be 
% relocated on a nearby voxel. The voxel to be relocated will be chosen 
% randomly from the neighbourhood, depending on how crowded it is (the 
% emptier a voxel is, the higher the probability to be chosen)

% Migration of a cell for a given genotype
%--------------------------------------------------------------------------
function gen = mig(gen,voxel,choice,Pmig)
    x = voxel(1);
    y = voxel(2);
    z = voxel(3);
    
    % We will be considering only adyacent movements
    movs = [-1,0,1];
    
    % Initially, cells will not migrate
    MIG = 0;
    
    r = rand;
    if r < Pmig
        MIG = 1;
    else
        MIG = 0;
    end
        
    if MIG == 1
        xmov = x+movs(randi(length(movs)));
        ymov = y+movs(randi(length(movs)));
        zmov = z+movs(randi(length(movs)));
        
        if xmov <= 1 || xmov >= size(gen,1)
            xmov = x;
        end
        
        if ymov <= 1 || ymov >= size(gen,1)
            ymov = y;
        end
        
        if zmov <= 1 || zmov >= size(gen,1)
            zmov = z;
        end

        gen{x,y,z}{choice,2} = gen{x,y,z}{choice,2}-1;
        
        % If a genotype is trying to move to a place where it does not
        % exists, it must be created before hand
        
        if size(gen{xmov,ymov,zmov},1) < choice
            nG = size(gen{xmov,ymov,zmov},1);
            gen{xmov,ymov,zmov}{choice,1} = gen{x,y,z}{choice,1};
            gen{xmov,ymov,zmov}{choice,2} = 1;
        else
            gen{xmov,ymov,zmov}{choice,2} = gen{xmov,ymov,zmov}{choice,2}+1;
        end
 
    end
    
    MIG = 0;
end

% Calculate probability of migration
%--------------------------------------------------------------------------
function Pmig = probmig(Hmig)
    Pmig = 1/(1+exp(-Hmig))-0.5;
end

% Calculate migration function
%--------------------------------------------------------------------------
function Hmig = migfit(voxel,gen,state,choice)
    % A 'voxel' contains all 3 coordinates needed to place it
    % 'gen' is the cell list containing all genotypes for each voxel
    % The 'state' is the state vector cell list for all voxels
    x = voxel(1);
    y = voxel(2);
    z = voxel(3);
    
    % Basic matrix for correlation between alterations (used to compute fitness) 
    M = zeros(6,6);
    M(1,1) = 0.25;
    M(2,2) = 0.25;
    M(3,3) = 0.25;
    
    nG = size(gen{x,y,z},1);
    genotype = gen{x,y,z}{choice,1};
    resources = state{x,y,z}(1)*state{x,y,z}(4);
    Hmig = genotype*M*genotype';
    for i = [1,nG]
        Gpop = gen{x,y,z}{i,2};
        Hmig = Hmig + 1/resources*Gpop;
    end

end


%% DEATH

% Everytime a death event takes place, the chosen genotype's population
% should be reduced in 1 individual

% Death of a given genotype
%--------------------------------------------------------------------------
function gen = death(gen,voxel,choice,Pkill)
    x = voxel(1);
    y = voxel(2);
    z = voxel(3);
    
    r = rand;
    if r < Pkill
        gen{x,y,z}{choice,2} = gen{x,y,z}{choice,2}-1;
    end
end

% Calculate probability of death
%--------------------------------------------------------------------------
function Pkill = probkill(voxel,gen,state,choice)
    x = voxel(1);
    y = voxel(2);
    z = voxel(3);
    
    Hkill = 0;
    nG = size(gen{x,y,z},1);
    genotype = gen{x,y,z}{choice,1};
    % 'resources' mixes both nutrients, angiogenesis, and limit population
    resources = state{x,y,z}(1)*state{x,y,z}(4);
    for i = [1,nG]
        Gpop = gen{x,y,z}{i,2};
        Hkill = Hkill + 1/resources*Gpop;
    end
    Pkill = 1/(1+exp(-Hkill))-0.5;
end


