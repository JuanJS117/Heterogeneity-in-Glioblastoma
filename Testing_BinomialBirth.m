P = 1e2;
K = 1e6;

T = [0:1:2*365*24/10];

Pop = zeros(1,length(T));
Pop(1) = P;
Preps = zeros(1,length(T));
Preps(1) = 0.01*(1-(Pop(t)/K));

for t = [1:length(T)-1]
    Prep = 0.01*(1-(Pop(t)/K));
    born = binornd(Pop(t),Prep);
    Pop(t+1) = Pop(t) + born;
    Preps(t+1) = Prep;
end

figure()
plot(T,Pop,'LineWidth',1)
xlabel('Time')
ylabel('Cell number')
title('Growth by binomial prob')

figure()
plot(T,Preps,'LineWidth',1)
xlabel('Time')
ylabel('Cell number')
title('Growth by binomial prob')