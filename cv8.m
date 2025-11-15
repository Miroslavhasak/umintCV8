clear; clc; close all;

u = (1.1*4+0.8*(4^3))/0.1; %  y = 4 y' = 0 y'' = 0
% y''+1.9*y'+1.1*y+0.8y^3=0.1*u

sim_name = ['cv8_sim'];

set_param(sim_name,'FastRestart','on');

tic
numgen=30	;% number of generations
lpop=30;	% number of chromosomes in population
lstring=6;	% number of genes in a chromosome       
M=10000;          % maximum of the search space


%parameters
a = 4; 
b = 2.3;
c = 0.0138;

Space=[zeros(1,lstring); ones(1,lstring)*M];

Delta=Space(2,:)/100;   

Pop=genrpop(lpop,Space);
%evolutions=zeros(1,numgen);

Fit = [];
Fit_a = [];
Fit_b = [];
Fit_c = [];

for gen=1:numgen

    for j=1:lpop
        disp(gen+"/"+numgen+" "+ j)
        P = Pop(j,1);
        I = Pop(j,2);
        D = Pop(j,3);

        P1 = Pop(j,4);
        I1 = Pop(j,5);
        D1 = Pop(j,6);

        try
            out=sim(sim_name);

            % pre prvy system
            Fit_1(j) = sum(a*abs(out.e.Data)+b*abs(out.de.Data)+c*abs(out.u.Data)); 

            % pre druhy system
            Fit_2(j) = sum(a*abs(out.e1.Data)+b*abs(out.de1.Data)+c*abs(out.u1.Data)); 

            % dokopy
            Fit(j) = Fit_1(j) + Fit_2(j);

        catch
            Fit(j) = 10e50;
        end
    end

%     GA
    Best=selbest(Pop,Fit,[1 1]);
    Old=selrand(Pop,Fit,9);
    Work1 = selsus(Pop,Fit,10);
    Work2 = selsus(Pop,Fit,9);
    Work1=crossov(Work1,1,0);
    Work2=mutx(Work2,0.1,Space);
    Work2=muta(Work2,0.2,Delta,Space);
    Pop=[Best;Old;Work1;Work2];

    evolution(gen)=min(Fit);
    evolution_1(gen) = min(Fit_1);
    evolution_2(gen) = min(Fit_2);

end
toc

plot(evolution,'k');
hold on
plot(evolution_1,'r');
plot(evolution_2,'g');
legend ('sum','1','2');
title('evolution');
xlabel('generation');
ylabel('fitness');
hold off


fprintf('Best parameters >> P: %.4f I: %.4f D: %.4f \n',Best(1,1),Best(1,2),Best(1,3));
fprintf('Best parameters >> P1: %.4f I1: %.4f D1: %.4f \n',Best(1,4),Best(1,5),Best(1,6));

P = Best(1,1);
I = Best(1,2);
D = Best(1,3);

P1 = Best(1,4);
I1 = Best(1,5);
D1 = Best(1,6);


figure
% pre hodnoty a
a0=1.1; a1=1.9; a2=1; b0=0.1; b1=0; c0=0.8;
out=sim(sim_name);
ya = out.y.Data;
plot(out.y.Data,'r');
hold on

% pre hodnoty b
a0=1.4; a1=2; a2=1; b0=0.2; b1=0; c0=0.8;
out=sim(sim_name);
yb = out.y.Data;
plot(out.y.Data,'g');


a0=1.1; a1=0.01; a2=4; b0=0.1; b1=0; c0=0.8;
out=sim(sim_name);
yc = out.y.Data;
plot(out.y.Data,'b');
plot(out.w.Data,'k');
title('priebeh regulovanej veliciny');
hold off

figure
% pre hodnoty a
a0=1.1; a1=1.9; a2=1; b0=0.1; b1=0; c0=0.8;
out=sim(sim_name);
ua = out.u.Data;
plot(ua,'r');
hold on

% pre hodnoty b
a0=1.4; a1=2; a2=1; b0=0.2; b1=0; c0=0.8;
out=sim(sim_name);
ub = out.u.Data;
plot(ub,'g');

% pre hodnoty c
a0=1.1; a1=0.01; a2=4; b0=0.1; b1=0; c0=0.8;
out=sim(sim_name);
uc = out.u.Data;
plot(uc,'b');
title('akcny zasah');

set_param(sim_name,'FastRestart','off');

% Best parameters >> P: 3481.9972 I: 1150.6154 D: 308.9246 
% Best parameters >> P1: 7564.2164 I1: 68.9133 D1: 628.8311 