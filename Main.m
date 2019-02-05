%% HETEROGENEITY IN GLIOBLASTOMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Preliminary program to test basic functionality of the model. It uses
% object-oriented programming to define a class (submodel) that represents
% a stores elemental information of the model in each lattice point.
%
% Fixed and genotype-independent reproduction, mutation, death and
% migration probabilities. 
%
% Associated files: submodel.m (class definition), reproduce.m, mutate.m,
% die.m
%
% -------------------------------------------------------------------------

tic
clear % clear workspace
clear all

% Define parameters

Xmax=15;
Ymax=15;

maxstep=100;

alt=2; % number of alterations
L=2^alt; % number of clones

% Define XxY lattice by creating Xmax * Ymax submodel objects (see definition).
% submodels is a Xmax by Ymax matrix 

submodels(1:Xmax,1:Ymax)=submod;

% Initialise submodels. Zero population, clones and existing clones.

for i=1:Xmax
    for j=1:Ymax
        
        submodels(i,j).N=0;
        submodels(i,j).clones=zeros(1,L);
        submodels(i,j).existent=[];
        
    end
end

% Place initial population (simplest initialisation). 10 individuals on
% each lattice point (or specificy initial point)

% for i=1:Xmax
%     for j=1:Ymax
i=8;j=8;
        submodels(i,j).clones(1) = 10; % place 10 individuals in genotype 1
        submodels(i,j).N=sum(submodels(i,j).clones); % compute N tot
        submodels(i,j).existent=[submodels(i,j).existent 1]; % add clone 1 
        % to existing clones
%     end
% end

evo=zeros(maxstep,2);
popin = sum([submodels.N]);
evo(1,2)=popin; % evo stores total population vs steps (to keep track of 
% variation in time)


area=zeros(maxstep,2);
area(1,2)=1; % area stores occupied lattice points vs steps (to keep track of 
% variation in time)

pop = zeros(Xmax,Ymax); % Array with total population on each point

clones = zeros(L,Xmax*Ymax); % Array with clonal population on each point

ex=zeros(L,Xmax*Ymax); % Array with existing clones on each point

% Begin update loop

submodelsnew=submodels; % auxiliar lattice

for step=1:maxstep % arbitrary number of steps
    
    % cycle through each lattice point for each step
    
    for i=1:Xmax
        for j=1:Ymax
            
            if submodels(i,j).N~=0 % operate only if lattice point is occupied
            
            chosen = datasample(submodels(i,j).existent,1); % Select individual at random
            
            % Reproduction %
            
            [submodelsnew] = reproduce(submodels,submodelsnew,i,j, chosen);
            
            % Mutation %
            
            [submodelsnew] = mutate(submodels,submodelsnew,i,j, chosen);
             
            % Death %
             
            [submodelsnew] = die(submodels,submodelsnew,i,j, chosen);
             
            % Migration 
             
            [submodelsnew] = migrate(submodels, submodelsnew,i,j, chosen, Xmax, Ymax);
             
            end
            
            pop(i,j)=submodels(i,j).N;
            
        end
    end
       
    submodels=submodelsnew; % Update
    
    % Compute total population
    
    poptot = sum([submodels.N]);
    
    % Save population and step in evo
    
    evo(step+1,:)=[step,poptot];
    area(step+1,:)=[step, nnz(pop)];
    
        % Plot after each 10% is completed
       
        if mod(step,maxstep/10)==0
    
            n=1;

        % extract information from objects and store in matrices

            for i=1:Xmax
                for j=1:Ymax
                    ex(1:length(submodels(i,j).existent),n)=submodels(i,j).existent;
                        clones(:,n)=submodels(i,j).clones;
                    pop(i,j)=submodels(i,j).N;
                     n=n+1;
                end
            end
        
            figure(step)

            imagesc(pop)
            title(['Population in space. Step', num2str(step)])
            set(gca, 'Fontsize',25,'LineWidth',3)
            colorbar 
            colormap(flipud(gray))
            caxis([0 ; max(max(pop))])
            print(['step ',num2str(step)],'-dpng','-r300')
        
        end
    
end


n=1;

% extract information from objects and store in matrices

for i=1:Xmax
    for j=1:Ymax
        ex(1:length(submodels(i,j).existent),n)=submodels(i,j).existent;
        clones(:,n)=submodels(i,j).clones;
        pop(i,j)=submodels(i,j).N;
        n=n+1;
    end
end

% output 

% pop
% clones
% ex

% Plot time evolution of population + linear fit

figure(1)

plot(evo(:,1),evo(:,2),'Linewidth',2)
set(gca, 'Fontsize',25,'LineWidth',3)
P = polyfit(evo(:,1),evo(:,2),1); % Linear fit
popfit = P(1)*evo(:,1)+P(2);
hold on;
plot(evo(:,1),popfit,'r-.','Linewidth',2);
title('Population vs time')
xlabel('Step','Fontsize',25)
ylabel('Population','Fontsize',25)
set(gca, 'Fontsize',25,'LineWidth',3)
hold off
print('evo','-dpng','-r300')
    
% Print slope and y-intercept

P(1)
P(2)

% Plot time evolution of volume

figure(2)

plot(area(:,1),area(:,2),'Linewidth',2)
set(gca, 'Fontsize',25,'LineWidth',3)
title('Area vs time')
xlabel('Step','Fontsize',25)
ylabel('Area','Fontsize',25)
print('area','-dpng','-r300')

toc


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------





% REPRODUCTION  -----------------------------------------------------------

function [ subup ] = reproduce(sub, subnew, n,m,chosen)

% Reproduction event with fixed probability
%
% Input
% -----
% sub = lattice to be updated
% subnew = transitory lattice 
% n,m = spatial coordinates of chosen individual
% chosen = genotype (integer) of the individual
%
% Output
% ------
% Updated population and clonal distribution after reproduction event
%
% -------------------------------------------------------------------------

if rand < 0.8
    
    subnew(n,m).clones(chosen)=sub(n,m).clones(chosen)+1; % add individual to clone
    subnew(n,m).N=sub(n,m).N+1; % increase total population
     
end
subup=subnew;
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% MUTATION  ---------------------------------------------------------------

function [  subup  ] = mutate(sub, subnew, n,m,chosen)

% Mutation event with fixed probability
%
% Input
% -----
% sub = lattice to be updated
% subnew = transitory lattice 
% n,m = spatial coordinates of chosen individual
% chosen = genotype (integer) of the individual
%
% Output
% ------
% Updated population and clonal distribution after mutation event
%
% -------------------------------------------------------------------------

if rand < 0.4
    
    subnew(n,m).clones(chosen)=subnew(n,m).clones(chosen)-1; % remove individual
    if subnew(n,m).clones(chosen)==0 % If no individuals remain, remove clone from 
        % existing clones
        subnew(n,m).existent = subnew(n,m).existent(subnew(n,m).existent~=chosen);
    end
    
    gen = de2bi(chosen-1,2); % Convert clone index to binary sequence 
    index = datasample(1:length(gen),1); % Choose element at random
    
    % Change element
    
    if gen(index)==1
        gen(index)=0;
    else
        gen(index)=1;
    end
    
    chosen = bi2de(gen)+1; % Convert back to integer
    
    subnew(n,m).clones(chosen)=subnew(n,m).clones(chosen)+1; % Add individual to mutated clone
    
    % If the mutated clone was unoccupied, add to list of existing clones
    
    if ismember(chosen,subnew(n,m).existent)==0 
        subnew(n,m).existent = [subnew(n,m).existent chosen];
    end
    
end
subup=subnew;
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% DEATH  ------------------------------------------------------------------

function [  subup  ] = die(sub, subnew, n,m,chosen)

% Death event with fixed probability
%
% Input
% -----
% sub = lattice to be updated
% subnew = transitory lattice 
% n,m = spatial coordinates of chosen individual
% chosen = genotype (integer) of the individual
%
% Output
% ------
% Updated population and clonal distribution after death event
%
% -------------------------------------------------------------------------

if subnew(n,m).clones(chosen)~=0 % Operate only if clone is still occupied after mutation

    if rand < 0.3
    
        subnew(n,m).clones(chosen)=subnew(n,m).clones(chosen)-1; % Remove individual from clone
        subnew(n,m).N=subnew(n,m).N-1; % Decrease total population
    
        % If no individuals remain, remove clone from existing clones
    
        if subnew(n,m).clones(chosen)==0
             subnew(n,m).existent = subnew(n,m).existent(subnew(n,m).existent~=chosen);
        end
     
    end

end
subup=subnew;
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% MIGRATION  --------------------------------------------------------------

function [subup] = migrate(sub, subnew, n, m, chosen, Xmax, Ymax)

% Migration event with fixed probability
%
% Input
% -----
% sub = lattice to be updated
% subnew = transitory lattice 
% n,m = spatial coordinates of chosen individual
% chosen = genotype (integer) of the individual
% Xmax, Ymax = lattice dimensions
%
% Output
% ------
% Updated population and clonal distribution after migration event
%
% -------------------------------------------------------------------------


if subnew(n,m).clones(chosen)~=0
    
    if rand < 0.8
    
         subnew(n,m).clones(chosen)=subnew(n,m).clones(chosen)-1; % remove individual
         subnew(n,m).N=subnew(n,m).N-1; % decrease total population
    
         % If no individuals left, remove clone from existing clones    
    
         if subnew(n,m).clones(chosen)==0
             subnew(n,m).existent = subnew(n,m).existent(subnew(n,m).existent~=chosen);
         end
    
         % Select migration destination (use info of previous lattice - sub)
    
         [i,j]=neighbour(sub,n,m,Xmax,Ymax);
    
         % Update population
    
         subnew(i,j).clones(chosen)=subnew(i,j).clones(chosen)+1; % add individual to clone
         subnew(i,j).N=subnew(i,j).N+1; % increase total population
    
         % If the mutated clone was unoccupied, add to list of existing clones
    
         if ismember(chosen,subnew(i,j).existent)==0 
            subnew(i,j).existent = [subnew(i,j).existent chosen];
         end
     
    end

end

subup=subnew;

end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% NEIGHBOUR SELECTION  ----------------------------------------------------

function [x,y] = neighbour(submodels,p,q,Xmax,Ymax)

% Select neighbour lattice point with weighted probability
%
% Input
% -----
% submodels = model in lattice point
% p,q = coordinates of chosen individual
% Xmax, Ymax = lattice dimensions
%
% Output
% ------
% x,y = coordinates of selected adjacent lattice point
%
%--------------------------------------------------------------------------

mat = zeros(8,3); % 8 possible destinations

index=0;

% cycle thru possible neighbours. Store coordinats if conditions are met
% (avoid boundaries)

for n=-1:1
    for m=-1:1
        
        l=p+n;
        k=q+m;
        
        if l<Xmax+1 && l>0 && k<Ymax+1 && k>0 && abs(n)+abs(m)~=0
            index=index+1;
            mat(index,:) = [l,k,submodels(l,k).N];
            
        end
        
    end
    
end

mat=mat(1:index,:); % Retain only accesible points

% Compute probabilities 

% mat(:,3)=1-mat(:,3)/sum(mat(:,3))

aux=sum(1./(mat(:,3)+1));

mat(:,3) = 1./(aux*(mat(:,3)+1));

% Select index in a with weighted probabilities

a=1:length(mat(:,1));

selection = a( sum( (rand(1) >= cumsum(mat(:,3)./sum(mat(:,3))))) + 1 );

 x=mat(selection,1);
 y=mat(selection,2);
            
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



