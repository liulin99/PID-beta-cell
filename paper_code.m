% to produce all Figures other than Figure 9 in the paper
names = ["G","G_S","G_I","G_Si","I","I_S","I_a","I_i","M_S","Ca_si","K_Si","Na_Si","V","A","sixCP"];
numbers = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
dic = dictionary(names,numbers);
options = odeset("NonNegative",[1     2     3     4     5     6     7     8     9    10    11    12    13    14    15]);
initial_glucose = 63.16; 

%the experimental point values are directly read from Figure 4G
time_points = [30,60,90,120,180,240,300,360];
T1D_glucose = [35,35,34, 34, 34, 32, 30.6,28.3];
Implant_glucose = [32.9,28,24.7,21.8,18.2,18.4,17,19.7];
Healthy = [14.5,9.74,8.24,8.6,7.56,6.26,5.43,4.14];


%% stage 1: 1 simulation
X0 = [0; 25; 0; 0; 0; 0; 0; 0; 0; 2*10^(-5); 200; 8;  -90; 0.1; 0]; %initial values used in stage 1
%simulate for three days without implant for healthy mice
tspan1 = 0:10:3*24*60;
[t, y1_1] = ode23(@(t,y) vivo_model(t,y,1,0,0,0), tspan1, X0); 
X1_1 = tail(y1_1,1);
% X1_1 = [5.94925519979674	25	5.94957447426638	0.0423894127678292	0.649186635621360	0	11.4656032357983	169.572292012184	0.0308740673470930	0.103499264959895	0.585284629098576	0.466763684964709	73.0680982513325	3.17207489695226	0.267904086138657];

%simulate for three days without implant for T1D mice
[t, y1_2] = ode23(@(t,y) vivo_model(t,y,0.001,0,0,0), tspan1, X0); 
X1_2 = tail(y1_2,1);
% X1_2 = [25.5861402490346	25	25.5862871761314	0.0423894127678292	0.00233550621044766	0	0.0411698696607124	94.7107645976048	0.0308740673470930	0.103499264959895	0.585284629098576	0.466763684964709	73.0680982513325	3.17207489695226	0.267904086138657];


%% stage 2: 2 simulations
tspan2 = 0:10:21*24*60;

X1_1(1,2) = 0; %set G_S to 0 
X1_2(1,2) = 0; %set G_S to 0

%simulation1: simulate for 21 days without implant for healthy mice
[t, y2_1] = ode23(@(t,y) vivo_model(t,y,1,0,0,0), tspan2, X1_1); 
X2_1 = tail(y2_1,1); %initial values for stage3 without implant
% X2_1 = [5.53723356890282	25	5.53723389333646	0.0423894127678293	0.710387799226937	0	12.5510163452660	157.756021217255	0.0308742396866240	0.103499263096041	0.585284413977824	0.466763624715294	73.0063888613527	3.17207489695228	0.267904086138662];

%simulation2: simulate for 21 days without implant for T1D mice
[t, y2_2] = ode23(@(t,y) vivo_model(t,y,0.001,0,0,0), tspan2, X1_2); 
X2_2 = tail(y2_2,1);
% X2_2 = [25.4483557332698	25	25.4483557437272	0.0423894127678293	0.00326526824898069	0	0.0576902450560541	507.126075856480	0.0308742396866240	0.103499263096041	0.585284413977824	0.466763624715294	73.0063888613527	3.17207489695228	0.267904086138662];

%simulation3: simulate for 21 days with a very low density(1 cell/mL)
%implant for sick mice to get the internal states of the beta cells for sick mice; then
%the blood insulin and glucose is substitute with the values from healthy fasting mice
[t, y2_3] = ode23(@(t,y) vivo_model(t,y,0.001,0.9137,0.7047,1), tspan2, X1_2); 
X2_3 = tail(y2_3,1);
% X2_3 = [25.6264216894093	25.6264216897280	25.6264216901862	0.0424158564766771	0.00189884114381980	0.00189655792181667	0.0335484296082557	510.257201386281	0.0308742395239717	0.103498772206185	0.586267230835652	0.466763031906142	73.0130618152451	3.17374557909211	0.268585211098882];

%simulation 4: simulate for 21 days with a density of 1.17*10^7 cell/mL,
%the density value is tuned to be visually consistent with Figure 4G
[t, y2_4] = ode23(@(t,y) vivo_model(t,y,0.001,0.9137,0.7047,1.17*10^7), tspan2, X1_2); 
X2_4 = tail(y2_4,1);
% X2_4 = [14.1644080480331	14.1644064298119	14.1644080633765	0.0415828657227478	0.159922813388294	0.160182642442112	2.82549137070896	304.656191178781	0.0308742443605029	0.103513356611723	0.557084573972080	0.466780887499238	73.0750934127015	3.12108831648045	0.247936152874373];

%% stage 3: 4 simulations
tspan3 = 0:1:6*60;
%set G_I to initial_glucose; then simulate for 6 hours of oral glucose
%tolerance test for 4 groups
X2_1(1,3) = initial_glucose;

X2_2(1,3) = initial_glucose;

X2_3(1,3) = initial_glucose;

X2_4(1,3) = initial_glucose;


%simulation1: healthy mice without any treatment(use X2_1 for initial values)
[t, y3_1] = ode23(@(t,y) vivo_model(t,y,1,0,0,0), tspan3, X2_1); 

%simulation2: T1D mice without any treatment(use X2_2 for initial values)
[t, y3_2] = ode23(@(t,y) vivo_model(t,y,0.001,0,0,0), tspan3, X2_2); 

%simulation3: T1D mice with high density beta cell treatment(use X2_3 for initial values)
[t, y3_3] = ode23(@(t,y) vivo_model(t, y, 0.001, 0.9137, 0.7047, 5*10^11), tspan3, X2_3);

%simulation4: T1D mice with high density PID beta cell treatment(use X2_3 for initial values)
[t, y3_4] = ode23(@(t,y) PID_model(t,y,0.001,0.9137,0.7047, 14, 11, 6, 3.2, 5*10^11), tspan3, [X2_3,0], options);

%simulation5: T1D mice with density of 1.17*10^7 for beta cell treatment(use X2_4 for initial values)
[t, low_density] = ode23(@(t,y) vivo_model(t,y,0.001,0.9137,0.7047,1.17*10^7), tspan3, X2_4); 

%use different target values
[t, y3_5] = ode23(@(t,y) PID_model(t,y,0.001,0.9137,0.7047, 14, 11, 6, 3.19, 5*10^11), tspan3, [X2_3,0], options);
[t, y3_6] = ode23(@(t,y) PID_model(t,y,0.001,0.9137,0.7047, 14, 11, 6, 3.21, 5*10^11), tspan3, [X2_3,0], options);

%use different PID gain parameters
[t, y3_7] = ode23(@(t,y) PID_model(t,y,0.001,0.9137,0.7047, 15, 12, 7, 3.2, 5*10^11), tspan3, [X2_3,0], options);
[t, y3_8] = ode23(@(t,y) PID_model(t,y,0.001,0.9137,0.7047, 13, 10, 5, 3.2, 5*10^11), tspan3, [X2_3,0], options);


%% result plotting

%reconstruct the Figure 4G
figure
hold on
plot(t,y3_1(:,dic("G")),'Color','k',LineWidth=2);
plot(t,y3_2(:,dic("G")),'Color','r',LineWidth=2);
plot(t,low_density(:,dic("G")),'Color','b',LineWidth=2);
scatter(time_points, Healthy,"black",'*');
scatter(time_points, T1D_glucose, "red",'o');
scatter(time_points, Implant_glucose,"blue",'filled');
legend({"WT","T1D","T1D A"});
xlabel('Time (minutes)');
ylabel('Glycemia (mM)');
hold off


%check the insulin and mRNA level in the cell
figure 
subplot(1,2,1)
hold on
plot(t,low_density(:,dic("I")),'Color','b',LineWidth=2);
plot(t,y3_1(:,dic("I")),'Color','k',LineWidth=2);
plot(t,y3_2(:,dic("I")),'Color','red',LineWidth=2);
% legend("T1D mice treated with artificial \beta-cells","healhty mice","T1D mice without treatment");
legend("T1D A","WT","T1D");
hold off
xlabel("Time (minutes)");
ylabel("Blood insulin (\mug/L)");
subplot(1,2,2)
hold on
plot(t,low_density(:,dic("M_S")),'Color','b',LineWidth=2);
xlabel("Time (minutes)");
ylabel("Intracellular mRNA concentration (mM)");
legend("T1D A");
ylim([0,0.05]);
hold off


figure
subplot(2,2,1)
hold on

plot(t,y3_1(:,dic("G")),'Color','k');
% plot(t,y3_2(:,dic("G")),'Color','r');
plot(t,y3_3(:,dic("G")),'Color','b','LineStyle',':',LineWidth=2);
% plot(t,y3_4(:,dic("G")),'Color','r',"LineStyle",'-.', LineWidth=2);

xlabel('Time (minutes)');
ylabel('Glycemia (mM)');

% legend({"WT","T1D","T1D A","T1D P"});
legend({"WT","T1D A"});

hold off


subplot(2,2,2)
hold on
plot(t,y3_3(:,dic("A")),'Color','b','LineStyle',':',LineWidth=2);
% plot(t,y3_4(:,dic("A")),'Color','r',"LineStyle",'-.',LineWidth=2);

xlabel('Time (minutes)');
ylabel('ATP concentration (mM)');
% legend({"T1D A","T1D P"});
legend({"T1D A"});
hold off

subplot(2,2,3)
hold on
plot(t,y3_1(:,dic("I")),'Color','k');
% plot(t,y3_2(:,dic("I")),'Color','r');
plot(t,y3_3(:,dic("I")),'Color','b','LineStyle',':',LineWidth=2); 
% plot(t,y3_4(:,dic("I")),'Color','r',"LineStyle",'-.',LineWidth=2);

xlabel('Time (minutes)');
ylabel('Blood Insulin (\mug/L)');
% legend({"WT","T1D","T1D A","T1D P"});
legend({"WT","T1D A"});
hold off

subplot(2,2,4)
hold on
% plot(t,y3_4(:,dic("M_S")),'Color','r',"LineStyle",'-.', LineWidth=2);
plot(t,y3_3(:,dic("M_S")),'Color','b','LineStyle',':',LineWidth=2);


xlabel('Time (minutes)');
ylabel('Cell mRNA concentration (mM)');
% legend({"T1D P","T1D A"});
legend({"T1D A"});
ylim([0,0.05]);
hold off

% disp(["distance compare",norm(y3_1(:,dic("G")) - y3_3(:,dic("G"))), norm(y3_1(:,dic("G")) - y3_4(:,dic("G")))]);



%for different PID gain parameter values
figure
hold on
plot(t, y3_4(:,dic("G")),'Color','m',LineWidth=2);
plot(t, y3_7(:,dic("G")),'Color','m',LineWidth=2,LineStyle=':');
plot(t, y3_8(:,dic("G")),'Color','m',LineWidth=2,LineStyle='-.');
legend({"P gain = 14, I gain = 11, D gain = 6","P gain = 15, I gain = 12, D gain = 7","P gain = 13, I gain = 10, D gain = 5"});
xlabel('Time (minutes)');
ylabel('Glycemia (mM)');
hold off


figure
hold on
plot(t,y3_4(:,dic("G")),'Color','m',LineWidth=2);
plot(t,y3_5(:,dic("G")),'Color','m',LineWidth=2,LineStyle='-.');
plot(t,y3_6(:,dic("G")),'Color','m',LineWidth=2,LineStyle=':');
xlabel('Time (minutes)');
ylabel('Glycemia (mM)');
legend({"ATP target value = 3.2","ATP target value = 3.19","ATP target value = 3.21"});
hold off



figure
subplot(2,2,1)
hold on

plot(t,y3_1(:,dic("G")),'Color','k');
plot(t,y3_3(:,dic("G")),'Color','b','LineStyle',':',LineWidth=2);
plot(t,y3_4(:,dic("G")),'Color','m',"LineStyle",'-.', LineWidth=2);

xlabel('Time (minutes)');
ylabel('Glycemia (mM)');

% legend({"WT","T1D","T1D A","T1D P"});
legend({"WT","T1D A","T1D P"});

hold off


subplot(2,2,2)
hold on
plot(t,y3_3(:,dic("A")),'Color','b','LineStyle',':',LineWidth=2);
plot(t,y3_4(:,dic("A")),'Color','m',"LineStyle",'-.',LineWidth=2);

xlabel('Time (minutes)');
ylabel('ATP concentration (mM)');
legend({"T1D A","T1D P"});
hold off
% disp(["distance compare",norm(y3_1(:,dic("G")) - y3_3(:,dic("G"))), norm(y3_1(:,dic("G")) - y3_4(:,dic("G")))]);

subplot(2,2,3)
hold on
plot(t,y3_1(:,dic("I")),'Color','k');
% plot(t,y3_2(:,dic("I")),'Color','r');
plot(t,y3_3(:,dic("I")),'Color','b','LineStyle',':',LineWidth=2); 
plot(t,y3_4(:,dic("I")),'Color','m',"LineStyle",'-.',LineWidth=2);

xlabel('Time (minutes)');
ylabel('Blood insulin (\mug/L)');
legend({"WT","T1D","T1D A","T1D P"});
legend({"WT","T1D A","T1D P"});
hold off


subplot(2,2,4)
hold on
plot(t,y3_4(:,dic("M_S")),'Color','m',"LineStyle",'-.', LineWidth=2);
plot(t,y3_3(:,dic("M_S")),'Color','b','LineStyle',':',LineWidth=2);


xlabel('Time (minutes)');
ylabel('Cell mRNA concentration (mM)');
legend({"T1D P","T1D A"});
hold off



%% ODE system definition
%gamma_G = 0.7047;  %/min, Vascular exchange constant glucose; Table S10
%gamma_I = 0.9137; %/min, Vascular exchange constant insulin; Table S10
%F_I: 0-1
function dVarsdt = vivo_model(t, Vars, F_I, gamma_I, gamma_G, N_S)
    %%constants
    g_p = 0.1061; %mM/min glucose production; Table S10
    g_r = 4.1*10^(-3); %mM/min glucose effectiveness, insulin-independent glucose removal; Table S10
    k_si = 1.2*10^(-3); %L/ug/min, Insulin sensitivity;  Table S10
    % V_S = 3.35*10^(-4); %L, Volume of the synthetic compartment for in-vivo experiments; Table S6
    % D_G = 0.2986; %/min, Vascular exchange constant glucose; Table S10
    % V_B = 1.6*10^(-3); %L, Volume of the mouse blood; Table S6
    % V_G = 3*10^(-3); %L, Volume of artificial glucose compartment; Table S6
    k_iar = 0.0566; %/min, Insulin clearance from action compartment; Table S10
    k_iir = 0.0351; %/min, Insulin clearance from virtual compartment; Table S10
    %F_I = 0.01; %factor for diabetic insulin production; Table S10 !!!
    % k_ir = 3.1439*10^(-4); %µg/L/mM/min, insulin clearance; Table S10
    V_max1 = 10^(-13); %mmol/min/cell, maximum specific rate of cellular glucose uptake; Table S7
    K_m1 = 0.31; %mM, Michaelis-Menten constant for the updtake of glucose; Table S7
    F_SCI = 0.0229; %ug/U, scaling factor of insulin; Table S6
    k_tl = 5*10^(-8); %U/mM/min/cell, common translation constant; Table S9
    k_di = 1.1*10^(-3); %/min, Degradation constant insulin in the synthetic compartment; Table S9
    n_NFAT = 9; %no unit, Number of NFAT repeats for calcium dependent Insulin expression; Table S9
    n_S = 3.6; %no unit, plasmid copy number for SEAP or Insulin transfection; Table S9 and S11
    n_C = 3.6; %no unit, plasmid copy number for Cav1.3 transfection; Table S9 and S11
    V_max4 = 2.67*10^(-6); %mM/min, Max. calcium-dependent gene expression; Table S9
    d_4 = 1.6; %no unit, Hill coefficient for calcium-dependent gene expression; Table S9, not found in the paper, already sent email for inqury; unclear !!!
    K_m4 = 1.1*10^(-3); %mM, affinity for calcium-dependent gene expression; Table S9
    k_dm = 2.8*10^(-3); %/min, Degradation constant mRNA; Table S9
    r = 6.5; %mu_m, Radius of HEK-293 cell; Table S6
    V_c = 3.1415926*4*(r*10^(-5))^3/3; %L, volume of HEK-293 cell; Table S6
    F = 9.6485*10^(4); %C/mol, Faraday constant; Table S8
    R = 8.3144; %J/mol/K, ideal gas constant; Table S8
    T = 310; %K, temperature; Table S8
    C = (7.2/60)*10^(-3); %nF, Capacitance of HEK-293 cell (scaled for time units in min.); Table S8; unclear !!!
    alpha = 6*10^(-8)/(V_c*F); %mM/min/pA, Coefficient converting current in pA to concentration per unit time in mM/min; Table S8
    f_Ca = 10^(-3); %no unit, Fraction of cytosolic calcium that is free; Table S8
    g_Ca = 6.43*10^(-3); %nS, Conductance of calcium ions through background calcium channels; Table S8
    g_CaV = 0.8; %nS, Conductance of calcium ions through voltage activated calcium channels; Table S8
    I_CaP_max = 9.4; %pA, Maximum calcium pump current; Table S8
    g_K = 1.97*10^(-2); %nS, conductance of potassium ions through background potassium channels; Table S8
    g_KV = 1.24; %nS, conductance of poatssium ions through voltage-activated potassium channels; Table S8
    g_KATP = 1.96; %nS, conductance of poatssium ions through ATP-dependent potassium channels; Table S8
    d_K = 7; %no unit, hill coefficient of ATP-sensitive potassium channles; Table S8
    A_tot = 5; %mM, total concentration of ATP+ADP; Table S7
    K1 = 3; %mM, ADP binding constant, ATP-sensitive potassium channels; Table S8
    K2 = 2.21; %mM, ATP binding constant, ATP-sensitive potassium channels; Table S8
    g_NA = 4.37*10^(-2); %nS, Conductance of sodium ions through background sodium channels; Table S8
    I_Nak_max = 42.72; %pA, Maximum sodiumpotassium pump curren; Table S8
    K_K = 1; %/min, Rate constant of potassium uptake from the medium by NaK pump; Table S8
    K_Na = 11; %/min, Rate constant of sodium uptake from the cell by NaK pump; Table S8
    q = 1; %no unit; Autocatalytic stoichiometry; Table S7
    u = 0.12; %mM/min; Max. global rate of ATP consumption; Table S7
    k_delta = 2; %mM; Michaelis-Menten constant for ATP consumption; Table S7
    K_m2A = 1.61; %mM; Affinity constant for allosteric inhibition of ATP consumption by ATP; Table S7
    a = 1;%ATP hill coefficient; Table S7
    b = 1;%ATP hill coefficient; Table S7
    c = 1;%ATP hill coefficient; Table S7
    V_max2 = 10; %mM/min; Max. rate of intracellular glucose and ATP consumption; Table S7
    K_m2G = 1.95; %mM; Michaelis-Menten constant for glucose consumption; Table S7
    V_max3 = 9.96; %mM/min; Max. intermediate consumption rate; Table S7
    K_m3 = 1.83; %mM; Affinity for intermediary consumption inhibition by ATP; Table S7
    Ca_S = 1.8; %mM; Table S4 and page 9, roughly constant
    K_S = 5.3; %mM; Table S4 and page 9, roughly constant
    Na_S = 155; %mM; Table S4 and page 9, roughly constant
    % k_p = 4.0338*10^(-5); %µg/L/mM/min; Basal insulin production; Table S10
    % k_i = 1.05*10^(-10); %µg/L/mM/min; Integral (insulin-mediated) insulin production; Table S10
    % k_d = 0.2340; %µg/L/mM/min; Differential (glucose-mediated) insulin production; Table S10
    V_S = 3.35*10^(-4); %L, Volume of the synthetic compartment for in-vivo experiments; Table S6
    %gamma_G = 0.7047;  %/min, Vascular exchange constant glucose; Table S10
    D_G = 0.2986; %/min, Vascular exchange constant glucose; Table S10
    V_B = 1.6*10^(-3); %L, Volume of the mouse blood; Table S6
    V_G = 3*10^(-3); %L, Volume of artificial glucose compartment; Table S6
    k_ir = 3.1439*10^(-4); %µg/L/mM/min, insulin clearance; Table S10
    k_p = 4.0338*10^(-5); %µg/L/mM/min; Basal insulin production; Table S10
    k_i = 1.05*10^(-10); %µg/L/mM/min; Integral (insulin-mediated) insulin production; Table S10
    k_d = 0.2340; %µg/L/mM/min; Differential (glucose-mediated) insulin production; Table S10
    %gamma_I = 0.9137; %/min; Vascular exchange constant insulin; Table S10
    %%end of constants
    dVarsdt = zeros(15,1);
    %define all 15 variables
    %Vars := [G,G_S,G_I,G_Si,I,I_S,I_a,I_i,M_S,Ca_si,K_Si,Na_Si,V,A,sixCP] ;
    G = Vars(1);
    G_S = Vars(2);
    G_I = Vars(3);
    G_Si = Vars(4);
    I = Vars(5);
    I_S = Vars(6);
    I_a = Vars(7);
    I_i = Vars(8);
    M_S = Vars(9);
    Ca_si = Vars(10);
    K_Si = Vars(11);
    Na_Si = Vars(12);
    V = Vars(13);
    A = Vars(14);
    sixCP = Vars(15);
    %copy all the necessary equations
    %[G,G_S,G_I,G_Si,I,I_S,I_a,I_i,M_S,Ca_si,K_Si,Na_Si,V,A,sixCP]
    %equation(1) is not needed in vivo; since mu will be zero in vivo(Table
    %S6),making mu_r(Ns) always zero in vivo
    % gamma_Gs = gamma_G*(G-G_S); %equation(37)
    gamma_GS = gamma_G*(G-G_S); %equation(37)
    R1_GS = G_S*V_max1/(G_S+K_m1); %equation(14)
    R_2_Gsi_A = (K_m2A^(2*b-a) * A^a * V_max2 * G_Si)/((K_m2A^(2*b)+A^(2*b))*(G_Si + K_m2G)); %Equation(15)
    dVarsdt(2) = -N_S*R1_GS + gamma_GS;%equation(2); dG_Sdt
    dVarsdt(4) = R1_GS/V_c - R_2_Gsi_A; %equation(3); dG_Sidt
    R_3_6CP_A = V_max3*(1-A/A_tot)*(K_m3^(2*c)*sixCP)/(A^(2*c+1)+K_m3^(2*c+1)); %Equation(16) !!!
    dVarsdt(15) = R_2_Gsi_A - R_3_6CP_A; %equation(4); dsixCPdt
    delta_A = u*A^2/(k_delta^2 + A^2); %Equation(18)
    dVarsdt(14) = -1*q*R_2_Gsi_A + (q+1)*R_3_6CP_A - delta_A; %equation(5); dAdt
    V_K_Ksi = 10^(3)*R*T*log(K_S/K_Si)/F; %equation(19)
    I_K_V_K_si = g_K*(V-V_K_Ksi); %equation(22)
    I_KV_V_K_si = g_KV*(V-V_K_Ksi) * ( 1/(1+exp((-15-V)/5.6)) ) * (1/(1+exp((43+V)/4.1))); %equation(23)
    I_KATP_V_KSi_A = g_KATP*(V-V_K_Ksi)*(1+ ((A_tot-A)/K1)^d_K )/(1+ ((A_tot-A)/K1)^d_K + (A/K2)^d_K); %equation(24)
    I_NaK_Na_Si = I_Nak_max*(K_S/(K_S+K_K))*(Na_Si/(Na_Si+K_Na))*((V+150)/(V+200)); %equation(28) unclear here without the V parameter, not defined in paper !!!
    I_NaK_V_Na_Si = I_Nak_max*(K_S/(K_S+K_K))*(Na_Si/(Na_Si+K_Na))*((V+150)/(V+200)); %equation(28)
    dVarsdt(11) = -1*alpha*(I_K_V_K_si + I_KV_V_K_si + I_KATP_V_KSi_A) + 2*alpha*I_NaK_Na_Si; %Equation(6) dK_Sidt
    V_Na_Nasi = 10^(3)*R*T*log(Na_S/Na_Si)/F; %equation(20)
    I_Na_V_NaSi = g_NA*(V-V_Na_Nasi); %equation(25)
    dVarsdt(12) = -1*alpha*I_Na_V_NaSi - 3*alpha*I_NaK_Na_Si; %Equation(7); dNa_sidt
    V_ca_Casi = 10^(3)*R*T*log(Ca_S/Ca_si)/F; %equation(21)
    I_Ca_V_Casi = g_Ca*(V-V_ca_Casi); %equation(26)
    I_Cav_V_Casi = n_C*g_CaV*(V - V_ca_Casi); %equation(27)
    I_CaP_V_Casi = I_CaP_max*Ca_si/(Ca_si+2*10^(-4));%equation(29)
    dVarsdt(10) = -0.5*alpha*f_Ca*(I_Ca_V_Casi+I_Cav_V_Casi+I_CaP_V_Casi); %Equation(8); dca_sidt
    dVarsdt(13) = (-1/C)*(I_K_V_K_si+I_KV_V_K_si+I_KATP_V_KSi_A+I_Na_V_NaSi+I_NaK_V_Na_Si+I_Ca_V_Casi+I_Cav_V_Casi+I_CaP_V_Casi); %Equation(9); dVdt
    R_4_CaSi = V_max4*(Ca_si)^d_4/(K_m4^d_4+(Ca_si)^d_4); %equation(17)
    dVarsdt(9) = n_NFAT*n_S*R_4_CaSi-(k_dm*M_S); %equation(10); dMsdt, mur(Ns) term is zero since mu is zero in vivo(Table S6)
    %equation(11); this is about the seap, not needed here
    gamma_IS = gamma_I*(I-I_S); %equation(36)
    dVarsdt(6) = F_SCI*k_tl*N_S*M_S - k_di*I_S + gamma_IS; %equation(12), dI_sdt
    %equation(13) is not needed since N_s will never change in vivo
    
    gamma_GI = D_G*(G-G_I); %equation(38)
    dVarsdt(1) = g_p - (g_r + k_si*I_a)*G - V_S*gamma_GS/V_B - V_B*gamma_GI/V_G; %equation(30), dGdt
    dVarsdt(3) = gamma_GI; %equation(31), dG_idt
    dVarsdt(7) = I - k_iar*I_a; %equation(32), dI_adt
    dVarsdt(8) = F_I*(G - k_iir*I_i); %equation(33), dI_idt, F_I value is unclear !!!
    n_HI = 8; %no unit
    k_MI = 10^(-2); %ng/ml in the original paper, which equals to µg/L in this paper
    cosi_I = I^n_HI/(I^n_HI + k_MI^n_HI); %unlabled, showed up on the very begining of page 14; unclear
    %Here is the guess:
    %https://www.cell.com/cms/10.1016/j.molcel.2014.06.007/attachment/992caabb-c4e2-4b4d-a3fb-9a0ddd8a1f60/mmc1.pdf
    %page 27
    %dVarsdt(5) = (F_I*(k_p*G + k_i*I_i + k_d*cosi_I*dVarsdt(1)) - V_S*gamma_IS/V_B - k_ir*I)*(1-I/20); %equation(34), dIdt; F_I value is unclear !!! add upper bount of 20mM for blood insulin
    dVarsdt(5) = (F_I*(k_p*G + k_i*I_i + k_d*cosi_I*dVarsdt(1)) - V_S*gamma_IS/V_B - k_ir*I)*(1-I/20);
    %equation(35), this is about the seap, not needed here
end

%%PID control with ATP
function dVarsdt = PID_model(t, Vars, F_I, gamma_I, gamma_G, p_gaining, i_gaining, d_gaining, steady_A, N_S)
    %%constants
    g_p = 0.1061; %mM/min glucose production; Table S10
    g_r = 4.1*10^(-3); %mM/min glucose effectiveness, insulin-independent glucose removal; Table S10
    k_si = 1.2*10^(-3); %L/ug/min, Insulin sensitivity;  Table S10
    V_S = 3.35*10^(-4); %L, Volume of the synthetic compartment for in-vivo experiments; Table S6
    D_G = 0.2986; %/min, Vascular exchange constant glucose; Table S10
    V_B = 1.6*10^(-3); %L, Volume of the mouse blood; Table S6
    V_G = 3*10^(-3); %L, Volume of artificial glucose compartment; Table S6
    k_iar = 0.0566; %/min, Insulin clearance from action compartment; Table S10
    k_iir = 0.0351; %/min, Insulin clearance from virtual compartment; Table S10
    %F_I = 0.01; %factor for diabetic insulin production; Table S10 !!!
    k_ir = 3.1439*10^(-4); %µg/L/mM/min, insulin clearance; Table S10
    V_max1 = 10^(-13); %mmol/min/cell, maximum specific rate of cellular glucose uptake; Table S7
    K_m1 = 0.31; %mM, Michaelis-Menten constant for the updtake of glucose; Table S7
    F_SCI = 0.0229; %ug/U, scaling factor of insulin; Table S6
    k_tl = 5*10^(-8); %U/mM/min/cell, common translation constant; Table S9
    k_di = 1.1*10^(-3); %/min, Degradation constant insulin in the synthetic compartment; Table S9
    n_NFAT = 9; %no unit, Number of NFAT repeats for calcium dependent Insulin expression; Table S9
    n_S = 3.6; %no unit, plasmid copy number for SEAP or Insulin transfection; Table S9 and S11
    n_C = 3.6; %no unit, plasmid copy number for Cav1.3 transfection; Table S9 and S11
    V_max4 = 2.67*10^(-6); %mM/min, Max. calcium-dependent gene expression; Table S9
    d_4 = 1.6; %no unit, Hill coefficient for calcium-dependent gene expression; Table S9, not found in the paper, already sent email for inqury; unclear !!!
    K_m4 = 1.1*10^(-3); %mM, affinity for calcium-dependent gene expression; Table S9
    k_dm = 2.8*10^(-3); %/min, Degradation constant mRNA; Table S9
    r = 6.5; %mu_m, Radius of HEK-293 cell; Table S6
    V_c = 3.1415926*4*(r*10^(-5))^3/3; %L, volume of HEK-293 cell; Table S6
    F = 9.6485*10^(4); %C/mol, Faraday constant; Table S8
    R = 8.3144; %J/mol/K, ideal gas constant; Table S8
    T = 310; %K, temperature; Table S8
    C = (7.2/60)*10^(-3); %nF, Capacitance of HEK-293 cell (scaled for time units in min.); Table S8; unclear !!!
    alpha = 6*10^(-8)/(V_c*F); %mM/min/pA, Coefficient converting current in pA to concentration per unit time in mM/min; Table S8
    f_Ca = 10^(-3); %no unit, Fraction of cytosolic calcium that is free; Table S8
    g_Ca = 6.43*10^(-3); %nS, Conductance of calcium ions through background calcium channels; Table S8
    g_CaV = 0.8; %nS, Conductance of calcium ions through voltage activated calcium channels; Table S8
    I_CaP_max = 9.4; %pA, Maximum calcium pump current; Table S8
    g_K = 1.97*10^(-2); %nS, conductance of potassium ions through background potassium channels; Table S8
    g_KV = 1.24; %nS, conductance of poatssium ions through voltage-activated potassium channels; Table S8
    g_KATP = 1.96; %nS, conductance of poatssium ions through ATP-dependent potassium channels; Table S8
    d_K = 7; %no unit, hill coefficient of ATP-sensitive potassium channles; Table S8
    A_tot = 5; %mM, total concentration of ATP+ADP; Table S7
    K1 = 3; %mM, ADP binding constant, ATP-sensitive potassium channels; Table S8
    K2 = 2.21; %mM, ATP binding constant, ATP-sensitive potassium channels; Table S8
    g_NA = 4.37*10^(-2); %nS, Conductance of sodium ions through background sodium channels; Table S8
    I_Nak_max = 42.72; %pA, Maximum sodiumpotassium pump curren; Table S8
    K_K = 1; %/min, Rate constant of potassium uptake from the medium by NaK pump; Table S8
    K_Na = 11; %/min, Rate constant of sodium uptake from the cell by NaK pump; Table S8
    q = 1; %no unit; Autocatalytic stoichiometry; Table S7
    u = 0.12; %mM/min; Max. global rate of ATP consumption; Table S7
    k_delta = 2; %mM; Michaelis-Menten constant for ATP consumption; Table S7
    K_m2A = 1.61; %mM; Affinity constant for allosteric inhibition of ATP consumption by ATP; Table S7
    a = 1;%ATP hill coefficient; Table S7
    b = 1;%ATP hill coefficient; Table S7
    c = 1;%ATP hill coefficient; Table S7
    V_max2 = 10; %mM/min; Max. rate of intracellular glucose and ATP consumption; Table S7
    K_m2G = 1.95; %mM; Michaelis-Menten constant for glucose consumption; Table S7
    V_max3 = 9.96; %mM/min; Max. intermediate consumption rate; Table S7
    K_m3 = 1.83; %mM; Affinity for intermediary consumption inhibition by ATP; Table S7
    Ca_S = 1.8; %mM; Table S4 and page 9, roughly constant
    K_S = 5.3; %mM; Table S4 and page 9, roughly constant
    Na_S = 155; %mM; Table S4 and page 9, roughly constant
    k_p = 4.0338*10^(-5); %µg/L/mM/min; Basal insulin production; Table S10
    k_i = 1.05*10^(-10); %µg/L/mM/min; Integral (insulin-mediated) insulin production; Table S10
    k_d = 0.2340; %µg/L/mM/min; Differential (glucose-mediated) insulin production; Table S10
    V_S = 3.35*10^(-4); %L, Volume of the synthetic compartment for in-vivo experiments; Table S6
    %gamma_G = 0.7047;  %/min, Vascular exchange constant glucose; Table S10
    D_G = 0.2986; %/min, Vascular exchange constant glucose; Table S10
    V_B = 1.6*10^(-3); %L, Volume of the mouse blood; Table S6
    V_G = 3*10^(-3); %L, Volume of artificial glucose compartment; Table S6
    k_ir = 3.1439*10^(-4); %µg/L/mM/min, insulin clearance; Table S10
    k_p = 4.0338*10^(-5); %µg/L/mM/min; Basal insulin production; Table S10
    k_i = 1.05*10^(-10); %µg/L/mM/min; Integral (insulin-mediated) insulin production; Table S10
    k_d = 0.2340; %µg/L/mM/min; Differential (glucose-mediated) insulin production; Table S10
    %gamma_I = 0.9137; %/min; Vascular exchange constant insulin; Table S10
    %%end of constants
    dVarsdt = zeros(16,1);
    %define all 15 variables
    %Vars := [G,G_S,G_I,G_Si,I,I_S,I_a,I_i,M_S,Ca_si,K_Si,Na_Si,V,A,sixCP] ;
    G = Vars(1);
    G_S = Vars(2);
    G_I = Vars(3);
    G_Si = Vars(4);
    I = Vars(5);
    I_S = Vars(6);
    I_a = Vars(7);
    I_i = Vars(8);
    M_S = Vars(9);
    Ca_si = Vars(10);
    K_Si = Vars(11);
    Na_Si = Vars(12);
    V = Vars(13);
    A = Vars(14);
    sixCP = Vars(15);

    %copy all the necessary equations
    %[G,G_S,G_I,G_Si,I,I_S,I_a,I_i,M_S,Ca_si,K_Si,Na_Si,V,A,sixCP]
    %equation(1) is not needed in vivo; since mu will be zero in vivo(Table
    %S6),making mu_r(Ns) always zero in vivo
    % gamma_Gs = gamma_G*(G-G_S); %equation(37)
    gamma_GS = gamma_G*(G-G_S); %equation(37)
    R1_GS = G_S*V_max1/(G_S+K_m1); %equation(14)
    R_2_Gsi_A = (K_m2A^(2*b-a) * A^a * V_max2 * G_Si)/((K_m2A^(2*b)+A^(2*b))*(G_Si + K_m2G)); %Equation(15)
    dVarsdt(2) = -N_S*R1_GS + gamma_GS;%equation(2); dG_Sdt
    dVarsdt(4) = R1_GS/V_c - R_2_Gsi_A; %equation(3); dG_Sidt
    %R_3_6CP_A = V_max3*(1-A/A_tot)*(K_m3^(2*c)*sixCP)/(A^(2*c)+K_m3^(2*c)); %Equation(16) !!!
    R_3_6CP_A = V_max3*(1-A/A_tot)*(K_m3^(2*c)*sixCP)/(A^(2*c+1)+K_m3^(2*c+1)); 
    dVarsdt(15) = R_2_Gsi_A - R_3_6CP_A; %equation(4); dsixCPdt
    delta_A = u*A^2/(k_delta^2 + A^2); %Equation(18)
    dVarsdt(14) = -1*q*R_2_Gsi_A + (q+1)*R_3_6CP_A - delta_A; %equation(5); dAdt
    V_K_Ksi = 10^(3)*R*T*log(K_S/K_Si)/F; %equation(19)
    I_K_V_K_si = g_K*(V-V_K_Ksi); %equation(22)
    I_KV_V_K_si = g_KV*(V-V_K_Ksi) * ( 1/(1+exp((-15-V)/5.6)) ) * (1/(1+exp((43+V)/4.1))); %equation(23)
    I_KATP_V_KSi_A = g_KATP*(V-V_K_Ksi)*(1+ ((A_tot-A)/K1)^d_K )/(1+ ((A_tot-A)/K1)^d_K + (A/K2)^d_K); %equation(24)
    I_NaK_Na_Si = I_Nak_max*(K_S/(K_S+K_K))*(Na_Si/(Na_Si+K_Na))*((V+150)/(V+200)); %equation(28) unclear here without the V parameter, not defined in paper !!!
    I_NaK_V_Na_Si = I_Nak_max*(K_S/(K_S+K_K))*(Na_Si/(Na_Si+K_Na))*((V+150)/(V+200)); %equation(28)
    dVarsdt(11) = -1*alpha*(I_K_V_K_si + I_KV_V_K_si + I_KATP_V_KSi_A) + 2*alpha*I_NaK_Na_Si; %Equation(6) dK_Sidt
    V_Na_Nasi = 10^(3)*R*T*log(Na_S/Na_Si)/F; %equation(20)
    I_Na_V_NaSi = g_NA*(V-V_Na_Nasi); %equation(25)
    dVarsdt(12) = -1*alpha*I_Na_V_NaSi - 3*alpha*I_NaK_Na_Si; %Equation(7); dNa_sidt
    V_ca_Casi = 10^(3)*R*T*log(Ca_S/Ca_si)/F; %equation(21)
    I_Ca_V_Casi = g_Ca*(V-V_ca_Casi); %equation(26)
    I_Cav_V_Casi = n_C*g_CaV*(V - V_ca_Casi); %equation(27)
    I_CaP_V_Casi = I_CaP_max*Ca_si/(Ca_si+2*10^(-4));%equation(29)
    dVarsdt(10) = -0.5*alpha*f_Ca*(I_Ca_V_Casi+I_Cav_V_Casi+I_CaP_V_Casi); %Equation(8); dca_sidt
    dVarsdt(13) = (-1/C)*(I_K_V_K_si+I_KV_V_K_si+I_KATP_V_KSi_A+I_Na_V_NaSi+I_NaK_V_Na_Si+I_Ca_V_Casi+I_Cav_V_Casi+I_CaP_V_Casi); %Equation(9); dVdt
    R_4_CaSi = V_max4*(Ca_si)^d_4/(K_m4^d_4+(Ca_si)^d_4); %equation(17)
    
    %dVarsdt(9) = n_NFAT*n_S*R_4_CaSi-(k_dm*M_S); %equation(10); dMsdt, mur(Ns) term is zero since mu is zero in vivo(Table S6)

    
    % steady_A = 3.08;%3.122852553601846; %get this value by 3days + 3 weeks + increase Gi and simulate for another 72 hours
    %tx1 = M_S-0.04;
    %cons = 10000;
    %dVarsdt(9) = 5*10^(-3)*(Vars(14)-steady_A) + 0*dVarsdt(14) + 0*(integral_error)-(k_dm*M_S)-(1+(exp(cons*tx1)-exp(-cons*tx1))/(exp(cons*tx1)+exp(-cons*tx1))); %mRNA is controlled by ATP  
    dVarsdt(16) = Vars(14)-steady_A; %the integration of ATP error
    %[Vars(9), Vars(14), (Vars(14)-steady_A), Vars(16), dVarsdt(14), (p_gaining*(Vars(14)-steady_A) + i_gaining*Vars(16) + d_gaining*(dVarsdt(14)) - (k_dm*M_S))*(1-Vars(9)/0.031)]
    dVarsdt(9) = (p_gaining*(Vars(14)-steady_A) + i_gaining*Vars(16) + d_gaining*(dVarsdt(14)) - (k_dm*M_S))*(1-Vars(9)/0.031);
    %dVarsdt(9) = p_gaining*(Vars(14)-steady_A) + i_gaining*Vars(16) + d_gaining*(dVarsdt(14)) - (k_dm*M_S);

    %equation(11); this is about the seap, not needed here
    gamma_IS = gamma_I*(I-I_S); %equation(36)
    dVarsdt(6) = F_SCI*k_tl*N_S*M_S - k_di*I_S + gamma_IS; %equation(12), dI_sdt
    %equation(13) is not needed since N_s will never change in vivo
    
    gamma_GI = D_G*(G-G_I); %equation(38)
    dVarsdt(1) = g_p - (g_r + k_si*I_a)*G - V_S*gamma_GS/V_B - V_B*gamma_GI/V_G; %equation(30), dGdt
    dVarsdt(3) = gamma_GI; %equation(31), dG_idt
    dVarsdt(7) = I - k_iar*I_a; %equation(32), dI_adt
    dVarsdt(8) = F_I*(G - k_iir*I_i); %equation(33), dI_idt, F_I value is unclear !!!
    n_HI = 8; %no unit
    k_MI = 10^(-2); %ng/ml in the original paper, which equals to µg/L in this paper
    cosi_I = I^n_HI/(I^n_HI + k_MI^n_HI); %unlabled, showed up on the very begining of page 14; unclear
    %Here is the guess:
    %https://www.cell.com/cms/10.1016/j.molcel.2014.06.007/attachment/992caabb-c4e2-4b4d-a3fb-9a0ddd8a1f60/mmc1.pdf
    %page 27
    tx2 = I-4.6383;
    %dVarsdt(5) = F_I*(k_p*G + k_i*I_i + k_d*cosi_I*dVarsdt(1)) - V_S*gamma_IS/V_B - k_ir*I - (1+(exp(cons*tx2)-exp(-cons*tx2))./(exp(cons*tx2)+exp(-cons*tx2)))*100 ; %equation(34), dIdt; F_I value is unclear !!!
    dVarsdt(5) = (F_I*(k_p*G + k_i*I_i + k_d*cosi_I*dVarsdt(1)) - V_S*gamma_IS/V_B - k_ir*I)*(1-I/20); %equation(34), dIdt; add upperbount of 20mM for blood insulin
    %equation(35), this is about the seap, not needed here
end
