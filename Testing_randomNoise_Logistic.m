%% TESTING PERTURBATION EFFECTS

close all; clear all; clc;

deltat = 10/24; K = 1e6; r = 1/48;

Ct = @(C,deltat,K,r) ((K*C*exp(r*deltat))/(K+C*(exp(r*deltat)-1)));

C0 = 1e2;

Nsteps = round(2*365*24/10);
T = linspace(0,Nsteps,Nsteps);

%% White noise in growth rate

C = zeros(1,Nsteps);
C(1) = C0;
Cpert1 = zeros(1,Nsteps);
Cpert1(1) = C0;
Cpert5 = zeros(1,Nsteps);
Cpert5(1) = C0;
Cpert10 = zeros(1,Nsteps);
Cpert10(1) = C0;

for t = [1:Nsteps-1]
    r1 = r*randn(1);
    r5 = r*(5+randn(1))/5;
    r10 = r*(10+randn(1))/10;
    C(t+1) = Ct(C(t),deltat,K,r);
    Cpert1(t+1) = Ct(Cpert1(t),deltat,K,r1);
    Cpert5(t+1) = Ct(Cpert5(t),deltat,K,r5);
    Cpert10(t+1) = Ct(Cpert10(t),deltat,K,r10);
end


figure()
hold on
plot(T,C,'LineWidth',1)
plot(T,Cpert1,'LineWidth',1)
plot(T,Cpert5,'LineWidth',1)
plot(T,Cpert10,'LineWidth',1)
hold off
legend('Without noise','Random pert. 1','Random pert. 5','Random pert. 10')


%% White noise in next step population

C = zeros(1,Nsteps);
C(1) = C0;
Cpert = zeros(1,Nsteps);
Cpert(1) = C0;


for t = [1:Nsteps-1]

    C(t+1) = Ct(C(t),deltat,K,r);
    Cpert(t+1) = Ct(Cpert(t),deltat,K,r);
    Cpert(t+1) = normrnd(Cpert(t+1),1e3*Cpert(t+1)/K);

end

figure()
hold on
plot(T,C,'LineWidth',1)
plot(T,Cpert,'LineWidth',1)
hold off
legend('Without noise','Random pert')





