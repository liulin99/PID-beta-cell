% to produce Figure 9 in the paper

density1 = 1.17*10^7;
density2 = 1.49*10^10;
density3 = 7.84*10^8;

names = ["G","G_S","G_I","G_Si","I","I_S","I_a","I_i","M_S","Ca_si","K_Si","Na_Si","V","A","sixCP"];
numbers = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
dic = dictionary(names,numbers);

time_points = [30,60,90,120,180,240,300,360];
T1D_glucose = [35,35,34, 34, 34, 32, 30.6,28.3];
Implant_glucose = [32.9,28,24.7,21.8,18.2,18.4,17,19.7];
Healthy = [14.5,9.74,8.24,8.6,7.56,6.26,5.43,4.14];

%%stage 1
X0 = [0; 25; 0; 0; 0; 0; 0; 0; 0; 2*10^(-5); 200; 8;  -90; 0.1; 0]; %initial values are copied from Table S4
%simulate for three days without implant for T1D mice
tspan = 0:10:3*24*60;
[t, y1] = ode23(@(t,y) vivo_model(t,y,0.001,0,0,1.17*10^7), tspan, X0); 
X1 = tail(y1,1);

%%stage 2
tspan = 0:10:21*24*60;
X1(1,2) = 0; %set G_S to 0
[t, y2_1] = ode23(@(t,y) vivo_model(t,y,0.001,0.9137,0.7047,density1), tspan, X1);
[t, y2_2] = ode23(@(t,y) vivo_model(t,y,0.001,0.9137,0.7047,density2), tspan, X1);
[t, y2_3] = ode23(@(t,y) vivo_model(t,y,0.001,0.9137,0.7047,density3), tspan, X1);
X2_1 = tail(y2_1,1);
X2_2 = tail(y2_2,1);
X2_3 = tail(y2_3,1);

%%stage 3
initial_glucose = 63.16; 
X2_1(1,3) = initial_glucose; 
X2_2(1,3) = initial_glucose; 
X2_3(1,3) = initial_glucose; 
tspan = 0:1:6*60;
[t, yd3_1] = ode23(@(t,y) vivo_model(t, y, 0.001, 0.9137, 0.7047, density1), tspan, X2_1);
[t, yd3_2] = ode23(@(t,y) vivo_model(t, y, 0.001, 0.9137, 0.7047, density2), tspan, X2_2);
[t, yd3_3] = ode23(@(t,y) vivo_model(t, y, 0.001, 0.9137, 0.7047, density3), tspan, X2_3);


figure
hold on
scatter(time_points, Implant_glucose,"blue",'filled');
plot(t,yd3_1(:,dic("G")),'Color','b',LineWidth=2);
plot(t,yd3_2(:,dic("G")),'Color','b',LineStyle=':',LineWidth=2);
plot(t,yd3_3(:,dic("G")),'Color','b',LineStyle='--',LineWidth=2);
legend({"Experiments mean blood glucose (8 mice)","Implanted cell density=1.17\times10^{7}","Implanted cell density=1.49\times10^{10}","Implanted cell density=7.84\times10^{8}"});
xlabel('Time (minutes)');
ylabel('Glycemia (mM)');
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
