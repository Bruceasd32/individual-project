%  State vector (14 states):
%    x(1) = ids   d-axis stator current       [A]
%    x(2) = iqs   q-axis stator current       [A]
%    x(3) = ldr   d-axis rotor flux           [Wb]
%    x(4) = lqr   q-axis rotor flux           [Wb]
%    x(5) = wr    rotor mechanical speed      [rad/s]
%    x(6) = id    d-axis grid current         [A]
%    x(7) = iq    q-axis grid current         [A]
%    x(8) = Vdc   dc-link voltage             [V]
%    x(9) = z1    generator-side ctrl state (= mds)
%    x(10)= z2    generator-side ctrl state (= mqs)
%    x(11)= z3    generator-side ctrl state (sphere closure)
%    x(12)= z4    grid-side ctrl state (= mdg)
%    x(13)= z5    grid-side ctrl state (= mqg)
%    x(14)= z6    grid-side ctrl state (sphere closure)
%
%  Controller is fully parameter-free; duty ratios are guaranteed to stay
%  inside the unit sphere (|mas|<=1, |mag|<=1) i.e. linear-modulation.
% -------------------------------------------------------------------------

clear; clc; close all;

%% ---- System parameters (Table I of the paper, 2 MW SCIG) ---------------
p  = struct();
p.Um   = 400*sqrt(2);      % peak phase voltage [V]  -> Vm
p.ws   = 100*pi;           % grid angular freq 50 Hz [rad/s]
p.Rg   = 0.05;             % boosting inductor resistance [Ohm]
p.Lg   = 1e-3;             % boosting inductor [H]
p.C    = 1e-3;             % dc-link capacitor [F]
p.Rdc  = 10e6;             % dc-link parasitic resistance [Ohm]
p.Rs   = 0.01;             % stator resistance [Ohm]
p.Ls   = 5.305e-3;         % stator inductance [H]
p.Rr   = 0.00842;          % rotor resistance [Ohm]
p.Lr   = 5.3137e-3;        % rotor inductance [H]
p.Lm   = 5.1839e-3;        % mutual inductance [H]
p.P    = 3;                % pole pairs
p.b    = 0.00015;          % friction [N.m.s/rad]
p.ng   = 62.5;             % gearbox ratio
p.Rb   = 35;               % blade radius [m]
p.Jwt  = 765.6;            % turbine inertia [kg.m^2]  (Jm in the paper)
p.rho  = 1.225;            % air density [kg/m^3]
p.sig  = p.Ls - p.Lm^2/p.Lr;   % leakage factor (sigma)

% Cp curve optimum (paper: lambda_opt ~ 6.325 for beta=0)
p.lambda_opt = 6.325;
p.beta       = 0;          % blade pitch angle

%% ---- References ---------------------------------------------------------
% Baseline values that have been verified to produce stable behaviour:
%   Vdc_ref = 1600 V, k1 = +0.03, k2 = +0.003, k3 = -1, k4 = -0.01
% (k1 and k3 signs are flipped relative to what is printed in the paper, see
%  the sign analysis in the original draft; Vdc_ref = 1600 V gives the duty
%  ratios enough headroom for the wind-step transients.)
p.ids_ref = 286;           % d-axis stator current ref [A]  (sets rotor flux)
p.Vdc_ref = 1600;          % dc-link voltage ref [V]

% Rotor time constant (used only inside the slip-freq formula)
p.tau_r = p.Lr / p.Rr;

%% ---- Controller gains  --------------------------------------------------
p.k1 = +0.03;                % [1/A]    ids   error gain
p.k2 = +0.003;               % [s/rad]  wr    error gain
p.k3 = -1;                   % [1/A]    id    error gain
p.k4 = -0.01;                % [1/V]    Vdc   error gain
p.c1 = 100;                  % sphere attractiveness (gen side)
p.c2 = 100;                  % sphere attractiveness (grid side)

%% ---- Wind-speed profile (Fig.3 of the paper) ---------------------------
%  11 m/s  ->  8 m/s (t=50) -> 10 m/s (t=150) -> 11.5 m/s (t=250)
p.wind_profile = @(t) ...
      11.0 *(t<20)  ...
    +  8.0 *(t>=20  & t<50) ...
    + 10.0 *(t>=50 & t<80) ...
    + 11.5 *(t>=80);

%% ---- Initial conditions: full analytical steady state at v=v0 ----------
% IMPORTANT: the closed-loop system is Lyapunov-stable, but the simulation
% is extremely stiff: at t=0 the back-emf (Lm/Lr)*p*lambda_dr*wr ~= 540 V
% must be balanced by 2*mqs*Vdc ~= 534 V, otherwise the q-axis stator
% voltage equation has a slope of order 1e6 A/s and the integrator drives
% the system out of the stability basin in milliseconds.  Therefore we
% solve the steady-state equations explicitly for iqs, iq and the four
% duty ratios so that all 14 derivatives are zero at t=0.
v0   = p.wind_profile(0);
wr0  = p.lambda_opt * v0 * p.ng / p.Rb;          % MPPT rotor speed
ldr0 = p.Lm * p.ids_ref;                          % from eq. (14)
lqr0 = 0;
Vdc0 = p.Vdc_ref;
ids0 = p.ids_ref;
id0  = 0;                                         % unity-PF objective

% --- mechanical equilibrium ---------------------------------------------
lam_i_opt = 1/(p.lambda_opt + 0.08*p.beta) - 0.035/(p.beta^3 + 1);
Cp_opt    = 0.22*(116*lam_i_opt - 0.4*p.beta - 5)*exp(-12.5*lam_i_opt);
Pwind0    = 0.5 * p.rho * pi * p.Rb^2 * v0^3;
Ptur0     = Cp_opt * Pwind0;
Tm0       = -Ptur0 / wr0;                         % wind drives the rotor
Te0       = p.b * wr0 + Tm0;                      % required electromag. torque

% --- iqs from Te = (3 Lm)/(2 Lr) * p * (ldr*iqs - lqr*ids), lqr = 0 ----
iqs0 = Te0 / ( 1.5 * (p.Lm/p.Lr) * p.P * ldr0 );

% --- generator-side duty ratios from steady-state stator equations ------
aRs0 = p.Rr*p.Lm^2/p.Lr^2 + p.Rs;
wsl0 = (1/(p.tau_r * p.ids_ref)) * iqs0;
we0  = p.P * wr0 + wsl0;
mds0 = ( aRs0*ids0 - p.sig*we0*iqs0 ...
        - (p.Rr*p.Lm/p.Lr^2)*ldr0 ) / (2*Vdc0);
mqs0 = ( aRs0*iqs0 + p.sig*we0*ids0 ...
        + (p.Lm/p.Lr)*p.P*ldr0*wr0 ) / (2*Vdc0);
z3_0 = sqrt(max(0, 1 - mds0^2 - mqs0^2));

% --- iq from dc-link power balance (Vdc_dot=0 with id=0)  --------------
%   3*mqg*iq = 3*(mds*ids + mqs*iqs) + Vdc/Rdc      and
%   2*mqg*Vdc = Vm - Rg*iq        =>  Rg*iq^2 - Vm*iq + 2*Vdc*S + 2*Vdc^2/(3 Rdc) = 0
S_gen = mds0*ids0 + mqs0*iqs0;
aq = p.Rg;  bq = -p.Um;  cq = 2*Vdc0*S_gen + 2*Vdc0^2/(3*p.Rdc);
iq0 = ( -bq - sqrt(bq^2 - 4*aq*cq) ) / (2*aq);    % delivering-power root
mqg0 = (p.Um - p.Rg*iq0) / (2*Vdc0);
mdg0 = (p.ws*p.Lg*iq0)   / (2*Vdc0);
z6_0 = sqrt(max(0, 1 - mdg0^2 - mqg0^2));

x0 = [ ids0;  iqs0;
       ldr0;  lqr0;
       wr0;
       id0;   iq0;
       Vdc0;
       mds0; mqs0; z3_0;
       mdg0; mqg0; z6_0 ];

%% ---- Simulate -----------------------------------------------------------
tspan = [0 120];
opts  = odeset('RelTol',1e-5,'AbsTol',1e-7,'MaxStep',0.005);

fprintf('Running simulation...\n');
tic;
[t,X] = ode23tb(@(t,x) scig_rhs(t,x,p), tspan, x0, opts);
fprintf('Done in %.1f s.  %d integration points.\n', toc, numel(t));

%% ---- Post-processing ----------------------------------------------------
ids = X(:,1);  iqs = X(:,2);
ldr = X(:,3);  lqr = X(:,4);
wr  = X(:,5);
id  = X(:,6);  iq  = X(:,7);
Vdc = X(:,8);
z1  = X(:,9);  z2  = X(:,10); z3 = X(:,11);
z4  = X(:,12); z5  = X(:,13); z6 = X(:,14);

% Modulation indices (eq. 4) -> must stay <= 1
mas = sqrt(z1.^2 + z2.^2);
mag = sqrt(z4.^2 + z5.^2);

% MPPT reference speed and wind-speed trace
v       = arrayfun(p.wind_profile, t);
wr_ref  = p.lambda_opt * v * p.ng / p.Rb;

% Powers injected into the grid (eq. 21/22 with the paper's sign convention
% Vd=0, Vq=Vm; Pg and Qg with -3/2 so that positive = generator->grid)
Vm  = p.Um;
Pg  = -1.5 * Vm .* iq;      % real power to grid (Vd=0)
Qg  = -1.5 * Vm .* id;      % reactive power to grid

% Mechanical / wind power
[Tm, Pwind] = arrayfun(@(vv,ww) turbine(vv,ww,p), v, wr);

%% ---- Recommended bounds from bounds_analysis.m -------------------------
% These four numbers come from running bounds_analysis.m (see the LaTeX
% document for full derivation).  They are drawn as black dashed lines on
% the four controlled-state subplots so that the simulation is visibly
% enveloped by the bounds.
beta_ids = 85.27;          % [A]      recommended bound for |e_ids|
beta_wr  = 33.88;          % [rad/s]  recommended bound for |e_wr|
beta_id  = 8.81;           % [A]      recommended bound for |e_id|
beta_Vdc = 94.06;          % [V]      recommended bound for |e_Vdc|

%% ---- Plots --------------------------------------------------------------
figure('Name','Wind speed','Position',[50 500 600 250]);
plot(t,v,'LineWidth',1.2); grid on;
xlabel('Time [s]'); ylabel('v [m/s]'); title('Wind speed profile');

figure('Name','SCIG states','Position',[50 50 1100 700]);
subplot(3,2,1); plot(t,ids,'LineWidth',1.1); hold on;
yline(p.ids_ref,'--r');
yline(p.ids_ref + beta_ids,'k--');
yline(p.ids_ref - beta_ids,'k--');
grid on; ylabel('i_{ds} [A]');
title(sprintf('(a) i_{ds} bound = ref \\pm %.2f A', beta_ids));
subplot(3,2,2); plot(t,iqs,'LineWidth',1.1); grid on;
ylabel('i_{qs} [A]'); title('(b) i_{qs} (free state)');
subplot(3,2,3); plot(t,ldr,'LineWidth',1.1); grid on;
ylabel('\lambda_{dr} [Wb]'); title('(c) \lambda_{dr}');
subplot(3,2,4); plot(t,lqr,'LineWidth',1.1); grid on;
ylabel('\lambda_{qr} [Wb]'); title('(d) \lambda_{qr} (near zero)');
subplot(3,2,5); plot(t,wr,'b', t,wr_ref,'--r','LineWidth',1.1); hold on;
plot(t, wr_ref + beta_wr,'k--', t, wr_ref - beta_wr,'k--');
grid on; ylabel('\omega_r [rad/s]'); xlabel('Time [s]');
legend('\omega_r','\omega_r^{ref}', sprintf('ref \\pm %.2f', beta_wr), ...
    'Location','best');
title('(e) \omega_r tracking with bound envelope');
subplot(3,2,6); plot(t,Vdc,'LineWidth',1.1); hold on;
yline(p.Vdc_ref,'--r');
yline(p.Vdc_ref + beta_Vdc,'k--');
yline(p.Vdc_ref - beta_Vdc,'k--');
grid on; ylabel('V_{dc} [V]'); xlabel('Time [s]');
title(sprintf('(f) V_{dc} bound = ref \\pm %.2f V', beta_Vdc));

figure('Name','Grid & modulation','Position',[700 500 900 500]);
subplot(2,2,1); plot(t,id,'LineWidth',1.1); hold on;
yline(+beta_id,'k--'); yline(-beta_id,'k--'); yline(0,'--r');
grid on; ylabel('i_d [A]');
title(sprintf('(g) i_d bound = \\pm %.2f A (unity PF)', beta_id));
subplot(2,2,2); plot(t,iq,'LineWidth',1.1); grid on;
ylabel('i_q [A]'); title('(h) q-axis grid current');
subplot(2,2,3); plot(t,mas,'LineWidth',1.1); grid on;
ylabel('m_{as}'); ylim([0 1.1]); yline(1,'--r');
title('(i) generator-side modulation index'); xlabel('Time [s]');
subplot(2,2,4); plot(t,mag,'LineWidth',1.1); grid on;
ylabel('m_{ag}'); ylim([0 1.1]); yline(1,'--r');
title('(j) grid-side modulation index'); xlabel('Time [s]');

figure('Name','Powers','Position',[700 50 900 400]);
subplot(1,2,1); plot(t,Pwind/1e6,'--k', t,Pg/1e6,'b','LineWidth',1.2);
grid on; xlabel('Time [s]'); ylabel('Power [MW]');
legend('P_{wind}','P_g','Location','best'); title('Active power');
subplot(1,2,2); plot(t,Qg/1e3,'LineWidth',1.2); grid on;
xlabel('Time [s]'); ylabel('Q_g [kVar]'); title('Reactive power (-> 0)');

%% ==========================================================================
function dx = scig_rhs(t, x, p)
% Right-hand-side of the 14-state closed-loop SCIG+converter+controller.
    ids = x(1); iqs = x(2); ldr = x(3); lqr = x(4); wr = x(5);
    id  = x(6); iq  = x(7); Vdc = x(8);
    z1  = x(9); z2  = x(10); z3 = x(11);
    z4  = x(12); z5 = x(13); z6 = x(14);

    % ---- References driven by wind --------------------------------------
    v      = p.wind_profile(t);
    wr_ref = p.lambda_opt * v * p.ng / p.Rb;        % eq. (9)

    % ---- Mechanical torque from the turbine Cp curve --------------------
    Tm = turbine(v, wr, p);

    % ---- Synchronous (ref-frame) speed via near-field-oriented slip -----
    % eq. (20): wsl = (1/(tau_r * ids_ref)) * iqs
    wsl = (1/(p.tau_r * p.ids_ref)) * iqs;
    we  = p.P * wr + wsl;                            % eq. (19)

    % ---- Controller duty ratios -----------------------------------------
    mds = z1;   mqs = z2;        % generator-side (eq. 15,16)
    mdg = z4;   mqg = z5;        % grid-side      (eq. 23,24)

    % ---- Induction generator (eq. 1) ------------------------------------
    sig  = p.sig;
    aRs  = p.Rr*p.Lm^2/p.Lr^2 + p.Rs;
    aLm  = p.Lm/p.Lr;

    dids = ( -aRs*ids + sig*we*iqs + (p.Rr*p.Lm/p.Lr^2)*ldr ...
             + aLm*p.P*lqr*wr + 2*mds*Vdc ) / sig;

    diqs = ( -aRs*iqs - sig*we*ids + (p.Rr*p.Lm/p.Lr^2)*lqr ...
             - aLm*p.P*ldr*wr + 2*mqs*Vdc ) / sig;

    dldr = (p.Rr*p.Lm/p.Lr)*ids - (p.Rr/p.Lr)*ldr + (we - p.P*wr)*lqr;
    dlqr = (p.Rr*p.Lm/p.Lr)*iqs - (p.Rr/p.Lr)*lqr - (we - p.P*wr)*ldr;

    % Mechanical equation (5th of eq.1): Jm*wr_dot = -1.5*(Lm/Lr)*p*lqr*ids
    %                                    + 1.5*(Lm/Lr)*p*ldr*iqs - b*wr - Tm
    Te   = 1.5*aLm*p.P*(ldr*iqs - lqr*ids);     % electromagnetic torque
    dwr  = (Te - p.b*wr - Tm) / p.Jwt;

    % ---- Grid-side converter + dc-link (eq. 2) --------------------------
    Vd = 0;  Vq = p.Um;                          % grid voltage alignment

    did  = (-p.Rg*id + p.ws*p.Lg*iq - 2*mdg*Vdc + Vd) / p.Lg;
    diq  = (-p.Rg*iq - p.ws*p.Lg*id - 2*mqg*Vdc + Vq) / p.Lg;
    dVdc = ( 3*(mdg*id + mqg*iq - mds*ids - mqs*iqs) - Vdc/p.Rdc ) / p.C;

    % ---- Nonlinear duty-ratio controllers (eq. 17,18 / 25,26) -----------
    e_ids = ids - p.ids_ref;
    e_wr  = wr  - wr_ref;
    s1    = z1^2 + z2^2 + z3^2 - 1;           % generator-side sphere residual

    dz1 = -p.k1 * e_ids * z3;
    dz2 = -p.k2 * e_wr  * z3;
    dz3 =  p.k1 * e_ids * z1 + p.k2 * e_wr * z2  -  p.c1 * s1 * z3;

    e_id  = id  - 0;                           % unity-PF -> id_ref = 0
    e_Vdc = Vdc - p.Vdc_ref;
    s2    = z4^2 + z5^2 + z6^2 - 1;

    dz4 = -p.k3 * e_id  * z6;
    dz5 = -p.k4 * e_Vdc * z6;
    dz6 =  p.k3 * e_id  * z4 + p.k4 * e_Vdc * z5  -  p.c2 * s2 * z6;

    dx = [dids; diqs; dldr; dlqr; dwr; did; diq; dVdc; ...
          dz1; dz2; dz3; dz4; dz5; dz6];
end

%% ==========================================================================
function [Tm, Pwind] = turbine(v, wr, p)
% Wind-turbine aerodynamic model (eq. 7,8).  Returns mechanical torque Tm
% reflected to the HIGH-SPEED (generator) shaft through the gearbox, and
% the ideal wind power Pwind = 0.5*rho*A*v^3 for reference.
    % tip-speed ratio referred to low-speed shaft
    v = max(v, 0.1);
    wr_ls = wr / p.ng;                      % turbine-shaft speed
    lambda = wr_ls * p.Rb / v;
    lam_i  = 1/(lambda + 0.08*p.beta) - 0.035/(p.beta^3 + 1);
    lami   = 1/lam_i;
    Cp     = 0.22*(116*lam_i - 0.4*p.beta - 5)*exp(-12.5*lam_i);
    Cp     = max(Cp, 0);                    % physical clamp

    A      = pi * p.Rb^2;
    Pwind  = 0.5 * p.rho * A * v^3;
    Ptur   = Cp * Pwind;                    % power captured by rotor
    % torque on the HIGH-speed (generator) shaft:
    wr_safe = max(abs(wr), 1e-3) * sign(wr + (wr==0));
    Tm      = -Ptur / wr_safe;
end

