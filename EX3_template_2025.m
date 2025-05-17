%% ME - 221 Dynamical systems, Spring 2025
% Solution to MATLAB PROBLEM SET 3
% Student Name:
% SCIPER:
% Submission date:

% Initiate the script
clc; clear; close all

%% Q6 - Compute the analytical parameters
% Define the parameters
m = 35*10^-3;  % kg
k = 16*10^5;  % N/m 
c = 32;  % Ns/m  

% Compute dimentionless parameters
w0 = sqrt(k/m);  % natural frequency
zeta = c/(2*sqrt(k*m));  % damping ratio
wbar = w0*sqrt(1-zeta^2);  % damped natural frequency

OS = 100*exp(-zeta*pi/sqrt(1-zeta^2));  % Maximum percent overshoot
tr = -sqrt(1-zeta^2)/(zeta*tan(wbar));  % rise time
ts = 4/(zeta*w0);  % 2% settling time

%% Q7a_tf Define the transfer function
sys_tf = tf(1/m, [1 c/m k/m]); %COMPLETE

%Compute the natural frequency and damping ratio using damp
[wn_tf,zeta_tf] = damp(sys_tf);

%Compute the overshoot, rise and settling time with stepinfo
TF_step_info = stepinfo(sys_tf,'SettlingTimeThreshold',0.02,'RiseTimeLimits',[0 1]);
OS_tf = TF_step_info.Overshoot;
tr_tf = TF_step_info.RiseTime;
ts_tf = TF_step_info.SettlingTime;

%% Q7a_ss Define the state-space model
% State is z = [y ydot]

A_ss = [0 1; -w0^2 -2*zeta*w0];
B_ss = [0; 1/m];
C_ss = [1 0];
D_ss = 0;

sys_ss = ss(A_ss, B_ss, C_ss, D_ss);

%Compute the natural frequency and damping ratio using damp
[wn_ss,zeta_ss] = damp(sys_ss);

%Compute the overshoot, rise and settling time with stepinfo
SS_step_info = stepinfo(sys_ss,'SettlingTimeThreshold',0.02,'RiseTimeLimits',[0 1]);
OS_ss = SS_step_info.Overshoot;
tr_ss = SS_step_info.RiseTime;
ts_ss = SS_step_info.SettlingTime;

%% Q7b Plot the step response
tt = linspace(0,1.2*ts_ss,1e4);

y_step = step(sys_ss, tt); %COMPLETE

COL = lines(4); COL(1,:) = [];
figure(1); clf
    plot(tt,y_step,'linewidth',1)
    hold all
    plot([SS_step_info.PeakTime SS_step_info.PeakTime],[0 SS_step_info.Peak],'--k','linewidth',0.5)
    plot([0 SS_step_info.PeakTime],[SS_step_info.Peak SS_step_info.Peak],'--k','linewidth',0.5)
    
    plot([0 tt(end)],1/k*[1 1],'--g','linewidth',0.5)
    plot([tr_ss tr_ss],1/k*[0 1],'--k','linewidth',0.5)
    
    plot([0 tt(end)],1/k*[1.02 1.02],'--k','linewidth',0.5)
    plot([0 tt(end)],1/k*[0.98 0.98],'--k','linewidth',0.5) 
    plot([ts_ss ts_ss],1/k*[0 1.02],'--k','linewidth',0.5) 
    
    h(1,1) = plot(SS_step_info.PeakTime,SS_step_info.Peak,'o','MarkerFaceColor',COL(1,:),'MarkerEdgeColor',COL(1,:));
    h(2,1) = plot(tr_ss,1/k,'o','MarkerFaceColor',COL(2,:),'MarkerEdgeColor',COL(2,:));
    h(3,1) = plot(ts_ss, 1/k*0.98,'o','MarkerFaceColor',COL(3,:),'MarkerEdgeColor',COL(3,:));
    
    hl = legend(h,'Over Shoot','Rise time','Settling time');
    set(hl,'interpreter','Latex')
    xlabel('$t$(s)','interpreter','Latex')
    ylabel('$y$(m)','interpreter','Latex')
    xlim([ 0 tt(end)])
    set(gca,'fontsize',8)
    set(gcf,'color','w','units','centimeters','position',[5 5 12 8])
    exportgraphics(gcf,'Figure_7.png','Resolution',300);

%% Q8 Material choice
E = [200 70 5]*1e9;   % Pa
rho = [7800 2700 1000]; % kg/m^3
alpha = [0.35 2 4];
beta = [5.9 1 8]*10^-4;
L = 0.2; % m
b = 0.1*L; % m
h = 0.2*L; % m
a = 0.9*L; % m

%% Q8a Plot the step response of the best combination of materials.
tt = linspace(0,3e-3,2e3).';

%Compute array with the parameters for the three materials
m = (b*h*a) * rho;
I = (1/12) * b * h^3;
ks = (3 * E * I) ./ (L^3 * (1 - a/L).^2);
c = alpha .* m + beta .* ks;
w0 = sqrt(ks ./ m); 
zetas = c ./ (2 .* sqrt(ks .* m));

tt = linspace(0,0.015,1e3);
for i=[1 2 3]
    sys = tf(1/m(i), [1 2*zetas(i)*w0(i) w0(i)^2]);
    step_info = stepinfo(sys,'SettlingTimeThreshold',0.02,'RiseTimeLimits',[0 1]);

    OS(i) = step_info.Overshoot;
    tr(i) = step_info.RiseTime;
    ts(i) = step_info.SettlingTime;
    Ss = step(sys,tt);
    y_step = step(sys,tt);
    
    COL = lines(4); COL(1,:) = [];
    figure; clf
    plot(tt,y_step,'linewidth',1)
    hold all
    plot([step_info.PeakTime step_info.PeakTime],[0 step_info.Peak],'--k','linewidth',0.5)
    plot([0 step_info.PeakTime],[step_info.Peak step_info.Peak],'--k','linewidth',0.5)
    
    plot([0 tt(end)],1/ks(i)*[1 1],'--k','linewidth',0.5)
    plot([tr(i) tr(i)],1/ks(i)*[0 1],'--k','linewidth',0.5)
    
    plot([0 tt(end)],1/ks(i)*[1.02 1.02],'--k','linewidth',0.5)
    plot([0 tt(end)],1/ks(i)*[0.98 0.98],'--k','linewidth',0.5) 
    plot([ts(i) ts(i)],1/ks(i)*[0 1.02],'--k','linewidth',0.5) 
    
    h(1,1) = plot(step_info.PeakTime,step_info.Peak,'o','MarkerFaceColor',COL(1,:),'MarkerEdgeColor',COL(1,:));
    h(2,1) = plot(tr(i),1/ks(i),'o','MarkerFaceColor',COL(2,:),'MarkerEdgeColor',COL(2,:));
    h(3,1) = plot(ts(i), 1/ks(i)*0.98,'o','MarkerFaceColor',COL(3,:),'MarkerEdgeColor',COL(3,:));
    
    hl = legend(h,'Over Shoot','Rise time','Settling time');
    set(hl,'interpreter','Latex','Location','SouthEast')
    xlabel('$t$(s)','interpreter','Latex')
    ylabel('$y$(m)','interpreter','Latex')
    xlim([ 0 tt(end)])
    set(gca,'fontsize',8)
    set(gcf,'color','w','units','centimeters','position',[5 5 12 8])
    exportgraphics(gcf,'Figure_8_'+string(i)+'.png','Resolution',300);
end

%% Q10 - Impulse response
M = 3; %kg
cd1 = 7; % N/sm
cd2 = 130; % N/sm
cd3 = 1650; % N/sm

A = 4500; % N
T = 10*10^-6; %s
f = 1/T;

% Selected material
id_mat = 2; % index of best material

m = m(id_mat)
k = ks(id_mat)
c = c(id_mat)

% Initial conditions
IC = zeros(4,1);
tt = linspace(0,2,1e5).';  % define the time vector
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

% Simulate the two cases           
[~,zNL_1P] = ode45(@(t,z)   EX3_NLode_2025(t,z,M,m,k,c,cd1,cd2,cd3,A,T), tt, IC, opts); % Positive impulse
[~,zNL_1N] = ode45(@(t,z)   EX3_NLode_2025(t,z,M,m,k,c,cd1,cd2,cd3,-A,T), tt, IC, opts); % Negative impulse

figure(6);clf
subplot(1,4,1:2)
    h = plot(tt,1e3*[zNL_1P(:,[1 2]) zNL_1N(:,[1 2])],'LineWidth',1.5);
    hold all
    plot([0 tt(end)],1e3*[zNL_1P(end,1) zNL_1P(end,1)],'--k','linewidth',0.5)
    plot([0 tt(end)],1e3*[zNL_1N(end,1) zNL_1N(end,1)],'--k','linewidth',0.5) 
    xlabel('$t$(s)','interpreter','Latex')
    ylabel('$x$(mm)','interpreter','Latex')
    hl = legend(h,'$y_P$','$u_P$','$y_N$','$u_N$');
    set(hl,'interpreter','Latex','Location','East')
    set(gca,'fontsize',8)
subplot(1,4,3:4)
    plot(tt,1e3*[zNL_1P(:,[1 2]) zNL_1N(:,[1 2])],'LineWidth',1.5);
    xlabel('$t$(s)','interpreter','Latex')
    xlim([0 0.6e-3])
    set(gca,'fontsize',8)
    set(gcf,'color','w','units','centimeters','position',[5 5 12 8])
    exportgraphics(gcf,'Figure_10.png','Resolution',300);

%% Q11 - Tracking performance

%Calculate the relative error
REP = (zNL_1P(:,2)-zNL_1P(:,1))./zNL_1P(:,2) * 100;
REN = (zNL_1N(:,2)-zNL_1N(:,1))./zNL_1N(:,2) * 100;
indP = find(abs(REP)<=10,1,'first');
indN = find(abs(REN)<=10,1,'first');

figure(7);clf
yyaxis left
    h = plot(tt,1e3*zNL_1P(:,[1 2]));
    xlabel('$t$(s)','interpreter','Latex')
    ylabel('$x$(mm)','interpreter','Latex')
    xlim([0 2e-3])
yyaxis right
    h(3,1) = plot(tt,REP); hold all
    plot(tt(indP)*[1 1],[-20 REP(indP)],'--k','LineWidth',1)
    plot([tt(indP) tt(end)],REP(indP)*[1 1],'--k','LineWidth',1)
    set(h(:),'LineWidth',1.5)
    hl = legend(h,'$y$','$u$','RE');
    set(hl,'interpreter','Latex','Location','East')
    ylabel('$RE$($\%$)','interpreter','Latex')
    xlim([0 2e-3])
    set(gca,'fontsize',8)
    set(gcf,'color','w','units','centimeters','position',[5 5 12 8])
    exportgraphics(gcf,'Figure_11a.png','Resolution',300);
figure(8);clf
yyaxis left
    h = plot(tt,1e3*zNL_1N(:,[1 2]));
    xlabel('$t$(s)','interpreter','Latex')
    ylabel('$x$(mm)','interpreter','Latex')
    xlim([0 2e-3])
yyaxis right
    h(3,1) = plot(tt,REN); hold all
    plot(tt(indN)*[1 1],[-20 REN(indN)],'--k','LineWidth',1)
    plot([tt(indN) tt(end)],REN(indN)*[1 1],'--k','LineWidth',1)
    set(h(:),'LineWidth',1.5)
    hl = legend(h,'$y$','$u$','RE');
    set(hl,'interpreter','Latex','Location','East')
    ylabel('$RE$($\%$)','interpreter','Latex')
    xlim([0 2e-3])
    set(gca,'fontsize',8)
    set(gcf,'color','w','units','centimeters','position',[5 5 12 8])
    exportgraphics(gcf,'Figure_11b.png','Resolution',300);
