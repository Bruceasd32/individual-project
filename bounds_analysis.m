% bounds_analysis.m
% Boundary range analysis of the four controlled state variables
%   ids, wr, id, Vdc
% in the closed-loop SCIG wind turbine system from
% scig_wt_nonlinear_control_2.m
%
% Computes four bounds and compares them:
%   B1 : empirical bound from full nonlinear ODE simulation
%   B2 : peak step response of the linearized closed-loop
%   B3 : Lyapunov ultimate bound
%   B4 : quasi-static reference jump + dynamic peak (initial-condition view)
%
% Outputs:
%   - console table
%   - bounds_analysis_results.txt
%   - bounds_analysis_results.mat
%   - several comparison figures

clear; clc; close all;

%% ---- Parameters --------------------------------------------------------
p = setupParams();

%% ---- Wind step set (paper Fig.3 profile) -------------------------------
windLevels = [11.0, 8.0, 10.0, 11.5];        % m/s
deltaV     = diff(windLevels);                % three steps
fprintf('Wind levels [m/s] : %s\n', mat2str(windLevels));
fprintf('Wind steps  [m/s] : %s\n', mat2str(deltaV));
fprintf('||Delta v||_inf   = %.2f m/s\n\n', max(abs(deltaV)));

%% ---- Equilibria at each wind level -------------------------------------
nLevels = numel(windLevels);
xEq = zeros(14, nLevels);
fprintf('=== Equilibrium points ===\n');
fprintf('%-6s %8s %8s %8s %8s %8s %8s %8s %8s\n', ...
    'v', 'iqs', 'wr', 'iq', 'mds', 'mqs', 'mdg', 'mqg', 'mas');
for j = 1:nLevels
    xEq(:, j) = computeEquilibrium(windLevels(j), p);
    mas = sqrt(xEq(9, j)^2 + xEq(10, j)^2);
    fprintf('%-6.2f %8.2f %8.3f %8.2f %8.4f %8.4f %8.4f %8.4f %8.4f\n', ...
        windLevels(j), xEq(2, j), xEq(5, j), xEq(7, j), ...
        xEq(9, j), xEq(10, j), xEq(12, j), xEq(13, j), mas);
end

%% ---- Linearization at each equilibrium ---------------------------------
% Closed-loop A and disturbance input matrix B (single input = wind speed v).
%
% Hurwitz evidence (addresses reviewer concern P2):
%   For each linearization point we report
%     - max real part of eig(A)        ->  must be < 0 for Hurwitz
%     - min damping ratio              ->  small value flags weakly damped mode
%     - cond(A)                        ->  large value flags ill conditioning
%     - ||m_s||^2, ||m_g||^2           ->  sphere-closure proximity to boundary
%   These data are dumped to A_metrics(:, :) and to the .txt results file.
A_all = cell(nLevels, 1);
B_all = cell(nLevels, 1);
A_metrics = zeros(nLevels, 5);   % cols: maxReEig, minDamp, condA, ||ms||^2, ||mg||^2
fprintf('\n=== Linearization evidence (Hurwitz + cond + sphere proximity) ===\n');
fprintf('%-6s %14s %12s %12s %12s %12s\n', ...
    'v', 'max Re(eig)', 'min zeta', 'cond(A)', '||m_s||^2', '||m_g||^2');
for j = 1:nLevels
    [A_all{j}, B_all{j}] = linearizeAtEq(xEq(:, j), windLevels(j), p);
    eigA   = eig(A_all{j});
    reEig  = real(eigA);
    imEig  = imag(eigA);
    absEig = abs(eigA);
    valid  = absEig > 1e-12;
    zetaAll = -reEig(valid) ./ absEig(valid);
    minZeta = min(zetaAll);
    condA   = cond(A_all{j});
    ms2 = xEq(9, j)^2  + xEq(10, j)^2;
    mg2 = xEq(12, j)^2 + xEq(13, j)^2;
    A_metrics(j, :) = [max(reEig), minZeta, condA, ms2, mg2];
    fprintf('%-6.2f %+14.3e %12.4f %12.3e %12.4f %12.4f\n', ...
        windLevels(j), max(reEig), minZeta, condA, ms2, mg2);
end
if all(A_metrics(:, 1) < 0)
    fprintf('  --> All linearization points are Hurwitz (max Re(eig) < 0).\n');
else
    warning('Some linearization point has max(Re(eig)) >= 0!');
end
fprintf('  --> Spheres far from boundary at all B2/B3 points: max ||m||^2 = %.4f << 1.\n', ...
    max(A_metrics(:, [4 5]), [], 'all'));

%% ---- Output matrices for the four controlled errors --------------------
% e = [e_ids; e_wr; e_id; e_Vdc] = C * delta_x + D * delta_v
% e_wr = (wr - wr_bar) - (lambda_opt*ng/Rb) * delta_v   -> direct feedthrough
C = zeros(4, 14);
C(1, 1) = 1;
C(2, 5) = 1;
C(3, 6) = 1;
C(4, 8) = 1;
D = [0; -p.lambda_opt * p.ng / p.Rb; 0; 0];

errLabels = {'|e_ids|  [A]', '|e_wr|   [rad/s]', '|e_id|   [A]', '|e_Vdc|  [V]'};
errKeys   = {'ids', 'wr', 'id', 'Vdc'};

%% ---- B1 : Empirical bound from full nonlinear simulation --------------
fprintf('\n=== B1 : Empirical bound from nonlinear simulation ===\n');
[tSim, xSim] = runNonlinearSim(p);
eSim = computeErrors(xSim, tSim, p);
beta_B1 = max(abs(eSim), [], 1)';
for k = 1:4
    fprintf('  %s : %.4g\n', errLabels{k}, beta_B1(k));
end

%% ---- B1 decomposition for wr channel (addresses reviewer concern P1) ----
% The wr-channel claim "B1 numerical truth equals B4 algebraic lower bound"
% is potentially tautological because both quantities are forced by the
% definition e_wr := wr - wr_ref(v) together with the unfiltered MPPT
% reference wr_ref(v) = 11.295*v: at any wind step Dv the reference jumps
% instantaneously by 11.295*Dv while wr (J_wt=765.6) cannot follow.
%
% We therefore separate the wr-channel peak into three independent parts:
%   (a) algebraic lower bound  alg_lb     = 11.295 * max|Dv|        (forced)
%   (b) instantaneous excursion right after each step
%       jump_at_step(j)        = |e_wr(t_step_j^+)|                 (forced)
%   (c) dynamic overshoot beyond the instantaneous excursion
%       dyn_overshoot(j)       = max_{t > t_step_j} |e_wr(t)|
%                                                   - jump_at_step(j)
% Claim "the controller dynamic does NOT add transient on wr" is supported
% only when dyn_overshoot is small relative to alg_lb. If dyn_overshoot is
% positive and non-negligible, the controller does add overshoot and the
% paper must state that explicitly.
fprintf('\n=== B1 decomposition for wr channel (P1 evidence) ===\n');
stepTimes_an = [20, 50, 80];
stepDvs      = diff(windLevels);
alg_lb_wr    = (p.lambda_opt * p.ng / p.Rb) * max(abs(stepDvs));

eWr = eSim(:, 2);
jumpVal = zeros(numel(stepTimes_an), 1);
ovrVal  = zeros(numel(stepTimes_an), 1);
fprintf('%-8s %14s %14s %14s\n', ...
    'step', '|e_wr(t+)|', 'peak after', 'overshoot');
for j = 1:numel(stepTimes_an)
    t0   = stepTimes_an(j);
    tNext = t0 + 30;   if j < numel(stepTimes_an), tNext = stepTimes_an(j + 1); end
    idxJust = find(tSim > t0, 1, 'first');
    idxWin  = find(tSim > t0 & tSim <= tNext);
    if isempty(idxJust) || isempty(idxWin)
        continue
    end
    jumpVal(j) = abs(eWr(idxJust));
    pkVal      = max(abs(eWr(idxWin)));
    ovrVal(j)  = max(0, pkVal - jumpVal(j));
    fprintf('Dv=%+4.1f  %14.4f %14.4f %14.4f\n', ...
        stepDvs(j), jumpVal(j), pkVal, ovrVal(j));
end
beta_wr_dynOvershoot = max(ovrVal);
beta_wr_jumpMax      = max(jumpVal);
fprintf('  algebraic lower bound  11.295 * max|Dv|         = %.4f rad/s\n', alg_lb_wr);
fprintf('  max instantaneous jump  |e_wr(t_step^+)|         = %.4f rad/s\n', beta_wr_jumpMax);
fprintf('  max controller-induced dynamic overshoot         = %.4f rad/s\n', beta_wr_dynOvershoot);
if beta_wr_dynOvershoot < 0.05 * alg_lb_wr
    fprintf('  --> Overshoot < 5%% of algebraic lower bound; controller\n');
    fprintf('      does not add transient beyond the definitional jump.\n');
else
    fprintf('  --> Overshoot is NON-negligible: must be reported separately.\n');
end

%% ---- B2 : Linearized step response peak --------------------------------
fprintf('\n=== B2 : Linearized step response peak ===\n');
beta_B2_per = zeros(4, numel(deltaV));
TendStep = 30;
for j = 1:numel(deltaV)
    A = A_all{j};   Bm = B_all{j};
    sysLin = ss(A, Bm, C, D);
    Ystep = step(deltaV(j) * sysLin, TendStep);
    beta_B2_per(:, j) = max(abs(Ystep), [], 1)';
    fprintf('  v=%4.1f -> %4.1f (Dv=%+4.1f) : %s\n', ...
        windLevels(j), windLevels(j+1), deltaV(j), ...
        sprintf('%10.4g', beta_B2_per(:, j)));
end
beta_B2 = max(beta_B2_per, [], 2);
fprintf('  ----- max over steps -----\n');
for k = 1:4
    fprintf('  %s : %.4g\n', errLabels{k}, beta_B2(k));
end

%% ---- B3a : Lyapunov ultimate bound (RAW, ill-scaled) -------------------
% Reproduces the pathological 3.8e15 result discussed in LaTeX section 5.2.
fprintf('\n=== B3a : Lyapunov ultimate bound, RAW (no scaling) ===\n');
W = max(abs(deltaV));
beta_B3raw = zeros(4, 1);
for j = 1:numel(windLevels)
    A = A_all{j};   Bm = B_all{j};
    Q = eye(14);
    P = lyap(A', Q);
    eP = eig(P);
    if min(eP) <= 0
        warning('P is not positive definite at j=%d; skipping.', j);
        continue
    end
    condFactor = sqrt(max(eP) / min(eP));
    PB_norm    = norm(P * Bm);
    lambdaQ    = min(eig(Q));
    ultBound   = condFactor * 2 * PB_norm * W / lambdaQ;
    for k = 1:4
        thisBound = norm(C(k, :)) * ultBound + abs(D(k)) * W;
        if thisBound > beta_B3raw(k)
            beta_B3raw(k) = thisBound;
        end
    end
end
for k = 1:4
    fprintf('  %s : %.4g  <-- pathological\n', errLabels{k}, beta_B3raw(k));
end

%% ---- B3b : Lyapunov ultimate bound, SCALED (LaTeX eq. T-29..T-31) ------
fprintf('\n=== B3b : Lyapunov ultimate bound, SCALED (Skogestad normalization) ===\n');
% Build the diagonal scaling S = diag(s_1,...,s_14) per LaTeX eq. (T-28).
% Each s_i = max( max_j |x_eq_i(v_j)|, nominal floor ).
floors = zeros(14, 1);
floors([1 2 6 7]) = 100;     % currents [A]
floors([3 4])      = 1;       % flux [Wb]
floors(5)          = 10;      % mechanical speed [rad/s]
floors(8)          = 100;     % dc-link voltage [V]
floors(9:14)       = 1;       % controller states z_i
S = max(max(abs(xEq), [], 2), floors);
fprintf('  Scaling factors s_i:\n');
fprintf('    %s\n', sprintf('%9.3g ', S));

beta_B3sc = zeros(4, 1);
errIdx = [1, 5, 6, 8];   % index of each controlled state in x
for j = 1:numel(windLevels)
    A     = A_all{j};
    Bm    = B_all{j};
    Sm    = diag(S);
    Sinv  = diag(1 ./ S);
    Atil  = Sinv * A * Sm;
    Btil  = Sinv * Bm;
    Q     = eye(14);
    Pt    = lyap(Atil', Q);
    eP    = eig(Pt);
    if min(eP) <= 0
        warning('Scaled P not p.d. at j=%d', j);
        continue
    end
    condFactor = sqrt(max(eP) / min(eP));
    PB_norm    = norm(Pt * Btil);
    ultTilde   = condFactor * 2 * PB_norm * W;       % ||delta_tilde_x||_ult
    for k = 1:4
        % per LaTeX eq. (T-31):  beta_*^{B3,sc} = s_{i_*} * ||delta_tilde_x|| + |d_*|*W
        thisBound = S(errIdx(k)) * ultTilde + abs(D(k)) * W;
        if thisBound > beta_B3sc(k)
            beta_B3sc(k) = thisBound;
        end
    end
    fprintf('  v=%4.1f : cond_factor=%9.3g, ||PB||=%9.3g, ||delta_tilde_x||_ult=%9.3g\n', ...
        windLevels(j), condFactor, PB_norm, ultTilde);
end
for k = 1:4
    fprintf('  %s : %.4g\n', errLabels{k}, beta_B3sc(k));
end

% Use the scaled B3 as the canonical "B3" entry in the comparison table.
beta_B3 = beta_B3sc;

%% ---- B4 : Quasi-static reference jump + dynamic peak -------------------
% After a step v_old -> v_new, the deviation from the NEW equilibrium has
% initial condition  delta_x(0+) = x_eq(:,old) - x_eq(:,new)
% The unforced linearized response gives the dynamic peak.
% For wr we additionally know the algebraic reference jump
%   |Delta wr_ref| = (lambda_opt*ng/Rb) * |Dv|
fprintf('\n=== B4 : Quasi-static + dynamic peak ===\n');
beta_B4 = zeros(4, 1);
TendIC = 30;
for j = 1:numel(deltaV)
    dx0   = xEq(:, j) - xEq(:, j + 1);
    A     = A_all{j + 1};
    sysIC = ss(A, zeros(14, 1), C, zeros(4, 1));
    Yic   = initial(sysIC, dx0, TendIC);
    peak  = max(abs(Yic), [], 1)';
    fprintf('  v=%4.1f -> %4.1f : %s\n', ...
        windLevels(j), windLevels(j + 1), sprintf('%10.4g', peak));
    for k = 1:4
        if peak(k) > beta_B4(k)
            beta_B4(k) = peak(k);
        end
    end
end
% Add the algebraic reference jump for wr:
algJump_wr = (p.lambda_opt * p.ng / p.Rb) * max(abs(deltaV));
beta_B4(2) = max(beta_B4(2), algJump_wr);
fprintf('  algebraic |Delta wr_ref| from max wind step = %.4f rad/s\n', algJump_wr);
fprintf('  ----- final B4 -----\n');
for k = 1:4
    fprintf('  %s : %.4g\n', errLabels{k}, beta_B4(k));
end

%% ---- Recommended bound -------------------------------------------------
% Take the maximum of empirical and analytical (B2/B4) so the bound is
% safe but not as conservative as B3.
betaRec = max([beta_B1, beta_B2, beta_B4], [], 2);

fprintf('\n=================================================================================\n');
fprintf('%-18s %10s %10s %10s %10s %10s %10s\n', ...
    'Error', 'B1', 'B2', 'B3raw', 'B3sc', 'B4', 'Recommended');
fprintf('%s\n', repmat('-', 1, 88));
for k = 1:4
    fprintf('%-18s %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n', ...
        errLabels{k}, beta_B1(k), beta_B2(k), beta_B3raw(k), beta_B3sc(k), ...
        beta_B4(k), betaRec(k));
end

%% ---- Absolute ranges ---------------------------------------------------
% References for ids, id, Vdc are constant; for wr the reference range is
% [min(wr_ref), max(wr_ref)] over the wind levels.
wrRefRange = [min(p.lambda_opt * p.ng / p.Rb * windLevels), ...
              max(p.lambda_opt * p.ng / p.Rb * windLevels)];
absRange = struct();
absRange.ids = [p.ids_ref - betaRec(1), p.ids_ref + betaRec(1)];
absRange.wr  = [wrRefRange(1) - betaRec(2), wrRefRange(2) + betaRec(2)];
absRange.id  = [-betaRec(3), betaRec(3)];
absRange.Vdc = [p.Vdc_ref - betaRec(4), p.Vdc_ref + betaRec(4)];

fprintf('\n=== Recommended absolute range of each controlled state ===\n');
fprintf('  ids in [%.2f, %.2f] A\n',     absRange.ids(1), absRange.ids(2));
fprintf('  wr  in [%.3f, %.3f] rad/s\n', absRange.wr(1),  absRange.wr(2));
fprintf('  id  in [%.4f, %.4f] A\n',     absRange.id(1),  absRange.id(2));
fprintf('  Vdc in [%.2f, %.2f] V\n',     absRange.Vdc(1), absRange.Vdc(2));

%% ---- Plots -------------------------------------------------------------
plotDir = pwd;
% Generate two sets of figures: Chinese (fig_*_zh.png) for bounds_theory.tex,
% English (fig_*_en.png) for bounds_theory_en.tex. Each set is 300 DPI.

% Bilingual title strings
titles_zh = struct( ...
    'bounds_linear',  '四种边界对比（线性坐标，已排除病态 B3 raw）', ...
    'bounds_log',     'Lyapunov 边界：病态 vs 归一化（对数坐标）', ...
    'sim_errors',     '14 阶非线性仿真误差与推荐边界包络', ...
    'theory_valid',   'B2 线性化阶跃预测 vs B1 非线性仿真（理论验证）', ...
    'wind_profile',   '风速廓线（论文 Fig.3）');
titles_en = struct( ...
    'bounds_linear',  'Comparison of bounds (linear scale, B3 raw excluded)', ...
    'bounds_log',     'Lyapunov bound: pathological vs scaled (log scale)', ...
    'sim_errors',     '14-state nonlinear simulation errors with recommended envelope', ...
    'theory_valid',   'B2 linearized prediction vs B1 nonlinear simulation', ...
    'wind_profile',   'Wind speed profile (Fig.3 of the paper)');
legend_zh = {'B1 仿真真值', 'B2 线性化预测', '推荐边界 \pm\beta'};
legend_en = {'B1 simulation truth', 'B2 linearized prediction', 'Recommended bound \pm\beta'};
recBoundLabel_zh = sprintf('推荐边界');
recBoundLabel_en = sprintf('recommended bound');

% --- Figure 1: linear-scale bar chart (excludes pathological B3raw) -----
makeBarLinear = @(ttl, fname) saveBarLinear(beta_B1, beta_B2, beta_B3sc, ...
    beta_B4, betaRec, errLabels, ttl, fullfile(plotDir, fname));
makeBarLinear(titles_zh.bounds_linear, 'fig_bounds_linear_zh.png');
makeBarLinear(titles_en.bounds_linear, 'fig_bounds_linear_en.png');

% --- Figure 2: log-scale comparison (shows pathological B3raw) ----------
makeBarLog = @(ttl, fname) saveBarLog(beta_B1, beta_B2, beta_B3raw, ...
    beta_B3sc, beta_B4, errLabels, ttl, fullfile(plotDir, fname));
makeBarLog(titles_zh.bounds_log, 'fig_bounds_log_zh.png');
makeBarLog(titles_en.bounds_log, 'fig_bounds_log_en.png');

% --- Figure 3: simulated errors with bound envelope ---------------------
makeSimEnv = @(ttl, lblBound, fname) saveSimEnv(tSim, eSim, betaRec, ...
    errLabels, errKeys, ttl, lblBound, fullfile(plotDir, fname));
makeSimEnv(titles_zh.sim_errors, recBoundLabel_zh, 'fig_sim_errors_zh.png');
makeSimEnv(titles_en.sim_errors, recBoundLabel_en, 'fig_sim_errors_en.png');

%% ---- Theory-vs-simulation validation ----------------------------------
% Reconstruct the predicted error trajectory by stitching together
%   * the linearized B2 step response after each wind step, and
%   * the algebraic jump for e_wr.
% Then plot it on top of the nonlinear B1 trajectory.
fprintf('\n=== Theory-vs-simulation validation ===\n');

stepTimes = [20, 50, 80];     % wind step instants
wrFactor  = p.lambda_opt * p.ng / p.Rb;

ePredict = zeros(size(eSim));   % zero before t=20 (initial equilibrium)
TpredictH = 30;                  % horizon for each linearized prediction

for j = 1:numel(stepTimes)
    A   = A_all{j};       % linearize at OLD wind level
    Bm  = B_all{j};
    sysLin = ss(A, Bm, C, D);
    tt = linspace(0, TpredictH, 1500)';
    Y  = step(deltaV(j) * sysLin, tt);     % 4-column output
    % Add this step response into ePredict at the appropriate times.
    timesAbs = stepTimes(j) + tt;
    mask     = timesAbs <= max(tSim);
    for k = 1:4
        ePredict(:, k) = ePredict(:, k) + ...
            interp1(timesAbs(mask), Y(mask, k), tSim, 'linear', 0) ...
            .* (tSim >= stepTimes(j));
        % subtract previous step's residual (only for currents/voltage that
        % return to zero anyway -- handled implicitly because step responses
        % converge before next step)
    end
end

% Compute prediction peaks and compare
beta_pred = max(abs(ePredict), [], 1)';
fprintf('%-18s %12s %12s %12s\n', 'Error', 'B1 truth', 'B2 predict', 'rel.err.');
for k = 1:4
    re = abs(beta_B1(k) - beta_pred(k)) / max(beta_B1(k), eps);
    fprintf('%-18s %12.4g %12.4g %12.2f%%\n', ...
        errLabels{k}, beta_B1(k), beta_pred(k), 100 * re);
end

% --- Figure 4: theory prediction overlaid on simulation ----------------
makeValid = @(ttl, leg, fname) saveValidation(tSim, eSim, ePredict, betaRec, ...
    errLabels, errKeys, beta_B1, beta_pred, ttl, leg, fullfile(plotDir, fname));
makeValid(titles_zh.theory_valid, legend_zh, 'fig_theory_validation_zh.png');
makeValid(titles_en.theory_valid, legend_en, 'fig_theory_validation_en.png');

% --- Figure 5: wind-speed profile --------------------------------------
makeWind = @(ttl, fname) saveWindProfile(tSim, p, ttl, fullfile(plotDir, fname));
makeWind(titles_zh.wind_profile, 'fig_wind_profile_zh.png');
makeWind(titles_en.wind_profile, 'fig_wind_profile_en.png');

fprintf('\nFigures saved (Chinese / English versions):\n');
for nm = {'wind_profile', 'bounds_linear', 'bounds_log', 'sim_errors', 'theory_validation'}
    fprintf('  fig_%s_zh.png   fig_%s_en.png\n', nm{1}, nm{1});
end

figure('Name', 'Empirical error trajectories', 'Position', [80 80 1100 700]);
for k = 1:4
    subplot(2, 2, k);
    plot(tSim, eSim(:, k), 'LineWidth', 1.0); hold on;
    yline(+betaRec(k), '--r', 'LineWidth', 1.0);
    yline(-betaRec(k), '--r', 'LineWidth', 1.0);
    grid on;
    xlabel('Time [s]');
    ylabel(errLabels{k});
    title(sprintf('%s   (recommended bound = %.4g)', errKeys{k}, betaRec(k)));
end
sgtitle('Nonlinear simulation errors vs. recommended bound');

%% ---- Save numerical results ----------------------------------=========
fid = fopen('bounds_analysis_results.txt', 'w');
fprintf(fid, '受控状态变量边界估计结果（自动生成）\n');
fprintf(fid, '================================================\n\n');
fprintf(fid, '风速廓线 [m/s] : %s\n', mat2str(windLevels));
fprintf(fid, '风速阶跃 [m/s] : %s\n', mat2str(deltaV));
fprintf(fid, '||Δv||_∞ = %.2f m/s\n\n', max(abs(deltaV)));

fprintf(fid, '--- 平衡点（v, iqs, wr, iq, mds, mqs, mdg, mqg, mas） ---\n');
for j = 1:nLevels
    mas = sqrt(xEq(9, j)^2 + xEq(10, j)^2);
    fprintf(fid, '%6.2f %10.3f %10.4f %10.3f %10.5f %10.5f %10.5f %10.5f %10.5f\n', ...
        windLevels(j), xEq(2, j), xEq(5, j), xEq(7, j), ...
        xEq(9, j), xEq(10, j), xEq(12, j), xEq(13, j), mas);
end

fprintf(fid, '\n--- 边界对比 ---\n');
fprintf(fid, '%-18s %10s %10s %10s %10s %10s %10s\n', ...
    'Error', 'B1', 'B2', 'B3raw', 'B3sc', 'B4', 'Recommended');
for k = 1:4
    fprintf(fid, '%-18s %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g\n', ...
        errLabels{k}, beta_B1(k), beta_B2(k), beta_B3raw(k), beta_B3sc(k), ...
        beta_B4(k), betaRec(k));
end

fprintf(fid, '\n--- 推荐绝对范围 ---\n');
fprintf(fid, '  ids ∈ [%.2f, %.2f] A\n',     absRange.ids(1), absRange.ids(2));
fprintf(fid, '  wr  ∈ [%.3f, %.3f] rad/s\n', absRange.wr(1),  absRange.wr(2));
fprintf(fid, '  id  ∈ [%.4f, %.4f] A\n',     absRange.id(1),  absRange.id(2));
fprintf(fid, '  Vdc ∈ [%.2f, %.2f] V\n',     absRange.Vdc(1), absRange.Vdc(2));

fprintf(fid, '\n--- 线性化点 Hurwitz / cond / 球面接近度 (P2 证据) ---\n');
fprintf(fid, '%-6s %14s %12s %12s %12s %12s\n', ...
    'v', 'max Re(eig)', 'min zeta', 'cond(A)', '||m_s||^2', '||m_g||^2');
for j = 1:nLevels
    fprintf(fid, '%-6.2f %+14.3e %12.4f %12.3e %12.4f %12.4f\n', ...
        windLevels(j), A_metrics(j, 1), A_metrics(j, 2), ...
        A_metrics(j, 3), A_metrics(j, 4), A_metrics(j, 5));
end

fprintf(fid, '\n--- wr 通道 B1 分解 (P1 证据) ---\n');
fprintf(fid, '  代数强制下界  11.295*max|Dv|             = %.4f rad/s\n', alg_lb_wr);
fprintf(fid, '  最大瞬时阶跃  max_j |e_wr(t_step_j^+)|    = %.4f rad/s\n', beta_wr_jumpMax);
fprintf(fid, '  最大控制器超调  max_j (peak - jump)        = %.4f rad/s\n', beta_wr_dynOvershoot);

fclose(fid);

save('bounds_analysis_results.mat', ...
    'windLevels', 'deltaV', 'xEq', 'A_all', 'B_all', 'A_metrics', ...
    'beta_B1', 'beta_B2', 'beta_B3raw', 'beta_B3sc', 'beta_B4', 'betaRec', ...
    'absRange', 'tSim', 'eSim', 'S', ...
    'alg_lb_wr', 'beta_wr_jumpMax', 'beta_wr_dynOvershoot', ...
    'jumpVal', 'ovrVal');

fprintf('\nResults saved to bounds_analysis_results.txt and .mat\n');

%% ========================================================================
%% Local helper functions
%% ========================================================================

function saveBarLinear(b1, b2, b3sc, b4, brec, lbl, ttl, fname)
% Bar chart, linear scale, 4 subplots (one per error)
    fig = figure('Position', [80 80 1100 700], 'Visible', 'off');
    for k = 1:4
        subplot(2, 2, k);
        bar([b1(k), b2(k), b3sc(k), b4(k), brec(k)]);
        set(gca, 'XTickLabel', {'B1', 'B2', 'B3sc', 'B4', 'Rec.'});
        title(lbl{k});
        grid on;
        ylabel('bound');
    end
    sgtitle(ttl);
    exportgraphics(fig, fname, 'Resolution', 300);
    close(fig);
end

function saveBarLog(b1, b2, b3raw, b3sc, b4, lbl, ttl, fname)
% Bar chart, log scale, 4 subplots, includes pathological B3raw
    fig = figure('Position', [80 80 1100 500], 'Visible', 'off');
    cats = categorical({'B1', 'B2', 'B3 raw', 'B3 scaled', 'B4'});
    cats = reordercats(cats, {'B1', 'B2', 'B3 raw', 'B3 scaled', 'B4'});
    for k = 1:4
        subplot(1, 4, k);
        bar(cats, [b1(k), b2(k), b3raw(k), b3sc(k), b4(k)]);
        set(gca, 'YScale', 'log');
        title(lbl{k});
        grid on;
        ylim([1e-2, 1e17]);
    end
    sgtitle(ttl);
    exportgraphics(fig, fname, 'Resolution', 300);
    close(fig);
end

function saveSimEnv(t, e, brec, lbl, keys, ttl, lblBound, fname)
% Time-domain error trajectories with bound envelope
    fig = figure('Position', [80 80 1100 700], 'Visible', 'off');
    for k = 1:4
        subplot(2, 2, k);
        plot(t, e(:, k), 'b', 'LineWidth', 1.0); hold on;
        yline(+brec(k), '--r', 'LineWidth', 1.0);
        yline(-brec(k), '--r', 'LineWidth', 1.0);
        grid on;
        xlabel('Time [s]');
        ylabel(lbl{k});
        title(sprintf('%s   (%s = %.4g)', keys{k}, lblBound, brec(k)));
    end
    sgtitle(ttl);
    exportgraphics(fig, fname, 'Resolution', 300);
    close(fig);
end

function saveValidation(t, e, ePred, brec, lbl, keys, b1, bp, ttl, leg, fname)
% Theory prediction overlaid on simulation
    fig = figure('Position', [80 80 1200 800], 'Visible', 'off');
    for k = 1:4
        subplot(2, 2, k);
        plot(t, e(:, k), 'b', 'LineWidth', 1.4); hold on;
        plot(t, ePred(:, k), 'r--', 'LineWidth', 1.0);
        yline(+brec(k), 'k:', 'LineWidth', 0.8);
        yline(-brec(k), 'k:', 'LineWidth', 0.8);
        grid on;
        xlabel('Time [s]');
        ylabel(lbl{k});
        legend(leg, 'Location', 'best', 'FontSize', 8);
        title(sprintf('%s : sim peak=%.3g, pred peak=%.3g', ...
            keys{k}, b1(k), bp(k)));
    end
    sgtitle(ttl);
    exportgraphics(fig, fname, 'Resolution', 300);
    close(fig);
end

function saveWindProfile(t, p, ttl, fname)
% Wind speed profile -- high-contrast styling (white background,
% filled step area, thick coloured step line, value labels).
    t = t(:);
    v = arrayfun(p.wind_profile, t);

    % crisp step coordinates from the piecewise-constant profile
    chg      = [true; diff(v) ~= 0];      % first sample of each plateau
    segStart = t(chg);
    segLevel = v(chg);
    segEnd   = [segStart(2:end); t(end)];
    tx = reshape([segStart, segEnd]',  [], 1);
    vy = reshape([segLevel, segLevel]', [], 1);

    lineCol = [0 0.30 0.63];      % deep blue step line
    fillCol = [0.80 0.90 1.00];   % light blue area
    mkCol   = [0.85 0.33 0.10];   % orange plateau markers
    yb      = 7.0;                % fill baseline [m/s]

    fig = figure('Position', [80 80 780 300], 'Visible', 'off', 'Color', 'w');
    ax  = axes(fig); hold(ax, 'on');

    fill([tx; flipud(tx)], [vy; yb*ones(size(vy))], fillCol, ...
         'EdgeColor', 'none', 'FaceAlpha', 0.9, 'Parent', ax);
    plot(ax, tx, vy, '-', 'Color', lineCol, 'LineWidth', 2.8);

    midT = (segStart + segEnd) / 2;
    plot(ax, midT, segLevel, 'o', 'MarkerSize', 7, ...
         'MarkerFaceColor', mkCol, 'MarkerEdgeColor', 'w', 'LineWidth', 1);
    for k = 1:numel(segLevel)
        text(ax, midT(k), segLevel(k) + 0.35, sprintf('%.1f', segLevel(k)), ...
            'Color', lineCol, 'FontWeight', 'bold', 'FontSize', 12, ...
            'HorizontalAlignment', 'center');
    end

    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'Color', 'w', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.6, ...
        'FontSize', 12, 'LineWidth', 1, 'Layer', 'top');
    xlabel(ax, 'Time [s]', 'FontSize', 13);
    ylabel(ax, 'v [m/s]',  'FontSize', 13);
    title(ax, ttl, 'FontSize', 13);
    xlim(ax, [t(1) t(end)]);
    ylim(ax, [yb, max(v) + 1]);

    exportgraphics(fig, fname, 'Resolution', 300);
    close(fig);
end

function p = setupParams()
% Reproduce the parameter struct from scig_wt_nonlinear_control_2.m
    p = struct();
    p.Um   = 400 * sqrt(2);
    p.ws   = 100 * pi;
    p.Rg   = 0.05;
    p.Lg   = 1e-3;
    p.C    = 1e-3;
    p.Rdc  = 10e6;
    p.Rs   = 0.01;
    p.Ls   = 5.305e-3;
    p.Rr   = 0.00842;
    p.Lr   = 5.3137e-3;
    p.Lm   = 5.1839e-3;
    p.P    = 3;
    p.b    = 0.00015;
    p.ng   = 62.5;
    p.Rb   = 35;
    p.Jwt  = 765.6;
    p.rho  = 1.225;
    p.sig  = p.Ls - p.Lm^2 / p.Lr;
    p.lambda_opt = 6.325;
    p.beta       = 0;
    p.ids_ref    = 286;
    p.Vdc_ref    = 1600;
    p.tau_r      = p.Lr / p.Rr;
    p.k1 = +0.03;
    p.k2 = +0.003;
    p.k3 = -1;
    p.k4 = -0.01;
    p.c1 = 100;
    p.c2 = 100;
    p.wind_profile = @(t) ...
          11.0 * (t < 20) ...
        +  8.0 * (t >= 20 & t < 50) ...
        + 10.0 * (t >= 50 & t < 80) ...
        + 11.5 * (t >= 80);
end

function x_eq = computeEquilibrium(v, p)
% Closed-form steady-state of the 14-state closed-loop system
% for a CONSTANT wind speed v.  Mirrors the analytical IC in the original
% script.
    wr   = p.lambda_opt * v * p.ng / p.Rb;
    ldr  = p.Lm * p.ids_ref;
    lqr  = 0;
    Vdc  = p.Vdc_ref;
    ids  = p.ids_ref;
    id   = 0;

    lam_i_opt = 1 / (p.lambda_opt + 0.08 * p.beta) - 0.035 / (p.beta^3 + 1);
    Cp_opt    = 0.22 * (116 * lam_i_opt - 0.4 * p.beta - 5) * exp(-12.5 * lam_i_opt);
    Pwind     = 0.5 * p.rho * pi * p.Rb^2 * v^3;
    Ptur      = Cp_opt * Pwind;
    Tm        = -Ptur / wr;
    Te        = p.b * wr + Tm;
    iqs       = Te / (1.5 * (p.Lm / p.Lr) * p.P * ldr);

    aRs  = p.Rr * p.Lm^2 / p.Lr^2 + p.Rs;
    wsl  = (1 / (p.tau_r * p.ids_ref)) * iqs;
    we   = p.P * wr + wsl;
    mds  = (aRs * ids - p.sig * we * iqs ...
            - (p.Rr * p.Lm / p.Lr^2) * ldr) / (2 * Vdc);
    mqs  = (aRs * iqs + p.sig * we * ids ...
            + (p.Lm / p.Lr) * p.P * ldr * wr) / (2 * Vdc);
    z3   = sqrt(max(0, 1 - mds^2 - mqs^2));

    S_gen = mds * ids + mqs * iqs;
    aq = p.Rg;
    bq = -p.Um;
    cq = 2 * Vdc * S_gen + 2 * Vdc^2 / (3 * p.Rdc);
    iq = (-bq - sqrt(bq^2 - 4 * aq * cq)) / (2 * aq);
    mqg = (p.Um - p.Rg * iq) / (2 * Vdc);
    mdg = (p.ws * p.Lg * iq) / (2 * Vdc);
    z6  = sqrt(max(0, 1 - mdg^2 - mqg^2));

    x_eq = [ids; iqs; ldr; lqr; wr; id; iq; Vdc; ...
            mds; mqs; z3; mdg; mqg; z6];
end

function [A, B] = linearizeAtEq(x_eq, v_eq, p)
% Numerical Jacobians of the closed-loop RHS at (x_eq, v_eq).
%   A = df/dx (14x14)    B = df/dv (14x1)
    n   = numel(x_eq);
    A   = zeros(n, n);
    eps = 1e-6;
    for i = 1:n
        xi = max(abs(x_eq(i)), 1.0) * eps;
        xp = x_eq;  xp(i) = xp(i) + xi;
        xm = x_eq;  xm(i) = xm(i) - xi;
        A(:, i) = (scigRhsXV(xp, v_eq, p) - scigRhsXV(xm, v_eq, p)) / (2 * xi);
    end
    dv = max(abs(v_eq), 1.0) * eps;
    B  = (scigRhsXV(x_eq, v_eq + dv, p) - scigRhsXV(x_eq, v_eq - dv, p)) / (2 * dv);
end

function dx = scigRhsXV(x, v, p)
% Same as scig_rhs in the original script but with wind speed v passed
% explicitly (no time dependence).
    ids = x(1);  iqs = x(2);  ldr = x(3);  lqr = x(4);  wr = x(5);
    id  = x(6);  iq  = x(7);  Vdc = x(8);
    z1  = x(9);  z2  = x(10); z3  = x(11);
    z4  = x(12); z5  = x(13); z6  = x(14);

    wr_ref = p.lambda_opt * v * p.ng / p.Rb;
    Tm     = turbineTorque(v, wr, p);

    wsl = (1 / (p.tau_r * p.ids_ref)) * iqs;
    we  = p.P * wr + wsl;

    mds = z1;  mqs = z2;
    mdg = z4;  mqg = z5;

    sig  = p.sig;
    aRs  = p.Rr * p.Lm^2 / p.Lr^2 + p.Rs;
    aLm  = p.Lm / p.Lr;

    dids = (-aRs * ids + sig * we * iqs + (p.Rr * p.Lm / p.Lr^2) * ldr ...
            + aLm * p.P * lqr * wr + 2 * mds * Vdc) / sig;
    diqs = (-aRs * iqs - sig * we * ids + (p.Rr * p.Lm / p.Lr^2) * lqr ...
            - aLm * p.P * ldr * wr + 2 * mqs * Vdc) / sig;
    dldr = (p.Rr * p.Lm / p.Lr) * ids - (p.Rr / p.Lr) * ldr ...
           + (we - p.P * wr) * lqr;
    dlqr = (p.Rr * p.Lm / p.Lr) * iqs - (p.Rr / p.Lr) * lqr ...
           - (we - p.P * wr) * ldr;
    Te   = 1.5 * aLm * p.P * (ldr * iqs - lqr * ids);
    dwr  = (Te - p.b * wr - Tm) / p.Jwt;

    Vd = 0;  Vq = p.Um;
    did  = (-p.Rg * id + p.ws * p.Lg * iq - 2 * mdg * Vdc + Vd) / p.Lg;
    diq  = (-p.Rg * iq - p.ws * p.Lg * id - 2 * mqg * Vdc + Vq) / p.Lg;
    dVdc = (3 * (mdg * id + mqg * iq - mds * ids - mqs * iqs) - Vdc / p.Rdc) / p.C;

    e_ids = ids - p.ids_ref;
    e_wr  = wr  - wr_ref;
    s1    = z1^2 + z2^2 + z3^2 - 1;
    dz1   = -p.k1 * e_ids * z3;
    dz2   = -p.k2 * e_wr  * z3;
    dz3   =  p.k1 * e_ids * z1 + p.k2 * e_wr * z2 - p.c1 * s1 * z3;

    e_id  = id;
    e_Vdc = Vdc - p.Vdc_ref;
    s2    = z4^2 + z5^2 + z6^2 - 1;
    dz4   = -p.k3 * e_id  * z6;
    dz5   = -p.k4 * e_Vdc * z6;
    dz6   =  p.k3 * e_id  * z4 + p.k4 * e_Vdc * z5 - p.c2 * s2 * z6;

    dx = [dids; diqs; dldr; dlqr; dwr; did; diq; dVdc; ...
          dz1; dz2; dz3; dz4; dz5; dz6];
end

function Tm = turbineTorque(v, wr, p)
% Wind-turbine aerodynamic torque on the high-speed shaft.
    v       = max(v, 0.1);
    wr_ls   = wr / p.ng;
    lambda  = wr_ls * p.Rb / v;
    lam_i   = 1 / (lambda + 0.08 * p.beta) - 0.035 / (p.beta^3 + 1);
    Cp      = 0.22 * (116 * lam_i - 0.4 * p.beta - 5) * exp(-12.5 * lam_i);
    Cp      = max(Cp, 0);
    A       = pi * p.Rb^2;
    Pwind   = 0.5 * p.rho * A * v^3;
    Ptur    = Cp * Pwind;
    wr_safe = max(abs(wr), 1e-3) * sign(wr + (wr == 0));
    Tm      = -Ptur / wr_safe;
end

function [t, X] = runNonlinearSim(p)
% Reproduce the full nonlinear simulation from the original script and
% return its solution.
    v0   = p.wind_profile(0);
    x0   = computeEquilibrium(v0, p);
    tspan = [0, 120];
    opts  = odeset('RelTol', 1e-5, 'AbsTol', 1e-7, 'MaxStep', 0.005);
    rhsFun = @(tt, x) scigRhsXV(x, p.wind_profile(tt), p);
    [t, X] = ode23tb(rhsFun, tspan, x0, opts);
end

function e = computeErrors(X, t, p)
% Return matrix of the four controlled errors at every sample time.
%   columns: e_ids, e_wr, e_id, e_Vdc
    ids = X(:, 1);
    wr  = X(:, 5);
    id  = X(:, 6);
    Vdc = X(:, 8);
    v       = arrayfun(p.wind_profile, t);
    wr_ref  = p.lambda_opt * p.ng / p.Rb * v;
    e = [ids - p.ids_ref, wr - wr_ref, id, Vdc - p.Vdc_ref];
end
