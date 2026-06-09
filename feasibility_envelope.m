%% feasibility_envelope.m
%  Concept A: Algebraic feasibility envelope of the SCIG closed-loop system.
%  No ODE integration is used to *find* the boundary; ODE is used only at
%  the very end as a one-shot verification at a boundary corner.
%
%  Methodology (matches the user's intuition):
%    1. List all physical / hardware constraints (Table 1 of the paper +
%       standard machine ratings):  v in [v_cut_in, v_cut_out], stator and
%       grid current ratings, DC-bus rectification floor, capacitor rating,
%       d-axis flux floor / saturation, mechanical speed limit.
%    2. List the controller-level constraints from the paper itself:
%       sphere constraints  m_ds^2 + m_qs^2 <= 1   and   m_dg^2 + m_qg^2 <= 1
%       (Lemma 1 of the original article; required for the duty-ratio
%       block to remain in the linear modulation region).
%    3. Sweep (v, Vdc_ref) over a wide rectangle. At each grid point compute
%       the closed-form equilibrium, evaluate every constraint, and mark the
%       cell feasible / which limit is binding.
%    4. Project the feasible set to the four controlled states ids, wr, id,
%       Vdc -> "state operating box".
%    5. Pick a corner of the feasibility boundary, plug the boundary
%       reference value into the original 14-state nonlinear model, run the
%       simulation, confirm the trajectory converges (does not diverge).
%
%  Output:  printed summary, .mat file, two PNG sets (Chinese/English).

close all; clear; clc;

p = setupParams();

%% ---- 1. Physical / hardware constraints -------------------------------
% Every value is derived from the original paper Table I + a published
% engineering standard / textbook.  See bounds_theory.tex Table 2 for
% the source citations corresponding to each row below.
LIM = struct();
LIM.v_min   = 4;                          % cut-in wind speed [m/s]
LIM.v_max   = 25;                         % cut-out wind speed [m/s]
% Vdc bounds: SVPWM modulation floor sqrt(6)*Vg and capacitor rating.
LIM.Vdc_min = sqrt(6) * 400;              % = 979.8 V (Yazdani 2010 Ch.6)
LIM.Vdc_max = 2400;                       % cap/IGBT rating, ~1.5 x nominal
% Stator current rating from machine nameplate calculation:
%   I_s,RMS = P_n / (sqrt(3) * V_LL * PF * eta)
%           = 2e6 / (sqrt(3) * 692.8 * 0.85 * 0.95) = 2063 A
%   I_s,max = sqrt(2) * 2063 = 2920 A (peak dq, amplitude-invariant Park)
% Sources: PF=0.85 from Bose, "Modern Power Electronics" Ch.2;
%          eta=0.95 from IEC 60034-2-1 (2-MW SCIG class).
LIM.is_max  = 2920;                       % stator winding rating [A peak]
% Grid converter rating from same nameplate at unity PF (PFC design):
%   I_g,RMS = P_n / (sqrt(3) * V_LL) = 2e6 / (sqrt(3) * 692.8) = 1667 A
%   I_g,max = sqrt(2) * 1667 = 2357 A peak.
LIM.ig_max  = 2357;                       % grid converter rating [A peak]
% ids floor / ceiling are absolute (not relative to ids_ref) because the
% closed-loop equilibrium gives ids = ids_ref.
% Floor: 0.3 x ids_ref = 86 A (avoid demagnetisation).
% Ceiling from rotor flux saturation, anchored on the *true* machine
% rated magnetising current i_m,rated = V_m / (omega_s * L_m):
%   i_m,rated = 565.69 / (100*pi * 5.1839e-3) = 347.4 A
%   lambda_rated = L_m * i_m,rated = V_m / omega_s = 1.801 Wb
%   lambda_sat   = 1.5 * lambda_rated = 2.702 Wb
%   i_ds,max,sat = lambda_sat / L_m = 2.702 / 5.1839e-3 = 521 A
% We adopt 425 A = 0.82 * 521 (~18 % safety margin under saturation).
% Reference: Krause/Wasynczuk/Sudhoff, Analysis of Electric Machinery, Ch.4.
% NOTE: this is *not* "1.5 * 286 = 429" -- the original paper labels
% i_ds_ref = 286 A as "rated", but it is actually a field-weakening
% operating point at 0.82 * i_m,rated.  Do not anchor on 286.
LIM.ids_min = 0.3 * 286;                  % = 86 A (anti-demagnetisation)
LIM.ids_max = 425;                        % = 0.82 * lambda_sat / L_m
% Mechanical speed limit on the high-speed shaft.  The full-power
% converter decouples generator electrical frequency from grid
% frequency, so the mechanical limit is set by the generator/gearbox
% nameplate, *not* by ws/P.  For a 2-MW class 6-pole SCIG, typical
% nameplate rated mechanical speed n_rated ~ 1500 rpm = 157 rad/s, with
% transient over-speed protection at 1.5-1.7 x n_rated (IEC 61400-1
% sec.11.4 trip threshold, ~235-267 rad/s).  We adopt 168 rad/s as a
% steady-state operating cap (~1.07 x n_rated, ~40-60% transient
% over-speed headroom).  The numerical coincidence with 1.6*ws/P is
% an artefact of common 6-pole 50-Hz nameplate practice, not a
% physical derivation.
LIM.wr_mech_max = 168;                    % steady-state mech cap [rad/s]
LIM.wr_min  = p.lambda_opt * LIM.v_min  * p.ng / p.Rb * 0.5;
LIM.wr_max  = LIM.wr_mech_max;

fprintf('=== Physical / hardware limits ===\n');
fprintf('  v       : [%.1f, %.1f] m/s\n',  LIM.v_min,   LIM.v_max);
fprintf('  Vdc_ref : [%.1f, %.1f] V\n',    LIM.Vdc_min, LIM.Vdc_max);
fprintf('  |i_s|   : <= %.0f A\n',         LIM.is_max);
fprintf('  |i_g|   : <= %.0f A\n',         LIM.ig_max);
fprintf('  ids     : [%.1f, %.1f] A\n',    LIM.ids_min, LIM.ids_max);
fprintf('  wr (MPPT) at v_max = %.1f rad/s\n', ...
    p.lambda_opt*LIM.v_max*p.ng/p.Rb);

%% ---- 2. Sweep (v, Vdc_ref) over a wide rectangle ---------------------
v_grid   = linspace(2, 30, 80);
Vdc_grid = linspace(400, 2800, 80);
nV = numel(v_grid); nD = numel(Vdc_grid);

feas      = false(nD, nV);
fail_code = zeros(nD, nV);
ms2_grid  = nan(nD, nV);
mg2_grid  = nan(nD, nV);
is_grid   = nan(nD, nV);
ig_grid   = nan(nD, nV);

% fail codes:
% 0 = feasible, 1 = wind range, 2 = Vdc range,
% 3 = stator sphere, 4 = grid sphere, 5 = |i_s|, 6 = |i_g|,
% 7 = wr_mech (rotor over-speed),
% -1 = numerical failure (e.g. complex sqrt)
fail_label = {'wind range', 'Vdc range', 'stator sphere ||m_s||>1', ...
              'grid sphere ||m_g||>1', '|i_s| > rating', '|i_g| > rating', ...
              'rotor over-speed wr > wr_mech'};

fprintf('\nSweeping %dx%d grid ...\n', nD, nV);
for iv = 1:nV
    for id_ = 1:nD
        v   = v_grid(iv);
        Vdc = Vdc_grid(id_);
        p_ij = p; p_ij.Vdc_ref = Vdc;

        if v < LIM.v_min || v > LIM.v_max
            fail_code(id_, iv) = 1; continue
        end
        if Vdc < LIM.Vdc_min || Vdc > LIM.Vdc_max
            fail_code(id_, iv) = 2; continue
        end

        try
            [eqv, ms2, mg2, ok_real] = computeEquilibriumDetail(v, p_ij);
        catch
            fail_code(id_, iv) = -1; continue
        end
        if ~ok_real
            fail_code(id_, iv) = -1; continue
        end

        ids_e = eqv(1); iqs_e = eqv(2); iq_e = eqv(7);
        is_mag = sqrt(ids_e^2 + iqs_e^2);
        ig_mag = abs(iq_e);
        ms2_grid(id_, iv) = ms2;
        mg2_grid(id_, iv) = mg2;
        is_grid(id_, iv)  = is_mag;
        ig_grid(id_, iv)  = ig_mag;

        wr_e = eqv(5);
        if ms2 > 1,                    fail_code(id_, iv) = 3; continue, end
        if mg2 > 1,                    fail_code(id_, iv) = 4; continue, end
        if is_mag > LIM.is_max,        fail_code(id_, iv) = 5; continue, end
        if ig_mag > LIM.ig_max,        fail_code(id_, iv) = 6; continue, end
        if wr_e   > LIM.wr_mech_max,   fail_code(id_, iv) = 7; continue, end

        feas(id_, iv) = true;
    end
end

%% ---- 3. Project to state operating boxes -----------------------------
[VV, DD] = meshgrid(v_grid, Vdc_grid);
idx      = find(feas);
v_feas   = VV(idx);
Vdc_feas = DD(idx);
wr_feas  = p.lambda_opt * v_feas * p.ng / p.Rb;

fprintf('\n=== Feasibility envelope summary ===\n');
fprintf('Sweep      : v in [%.1f, %.1f] m/s, Vdc_ref in [%.0f, %.0f] V\n', ...
    v_grid(1), v_grid(end), Vdc_grid(1), Vdc_grid(end));
fprintf('Feasible cells : %d / %d  (%.1f%%)\n', ...
    numel(idx), nV*nD, 100*numel(idx)/(nV*nD));
fprintf('  v       feasible in [%.2f, %.2f] m/s\n', min(v_feas), max(v_feas));
fprintf('  Vdc_ref feasible in [%.1f, %.1f] V\n',   min(Vdc_feas), max(Vdc_feas));
fprintf('  wr_eq   feasible in [%.2f, %.2f] rad/s (= MPPT(v))\n', ...
    min(wr_feas), max(wr_feas));

% Round to grid resolution (the sweep's discretisation step) for honest
% reporting; the binding constraint (e.g. is_max, wr_mech_max) gives the
% physical edge, the grid cell is just the resolution we sampled at.
dv  = v_grid(2)   - v_grid(1);
dDV = Vdc_grid(2) - Vdc_grid(1);

BOX = struct();
BOX.ids_lo = p.ids_ref;       BOX.ids_hi = p.ids_ref;
BOX.wr_lo  = min(wr_feas);    BOX.wr_hi  = max(wr_feas);
BOX.id_lo  = 0;               BOX.id_hi  = 0;
BOX.Vdc_lo = min(Vdc_feas);   BOX.Vdc_hi = max(Vdc_feas);
BOX.v_lo   = min(v_feas);     BOX.v_hi   = max(v_feas);

fprintf('\n=== State operating boxes (Concept A, steady-state envelope) ===\n');
fprintf('  v   in [%7.2f, %7.2f]  m/s    (sweep step %.2f m/s)\n', ...
    BOX.v_lo, BOX.v_hi, dv);
fprintf('  ids in [%7.2f, %7.2f]  A      (fixed reference)\n', ...
    BOX.ids_lo, BOX.ids_hi);
fprintf('  wr  in [%7.2f, %7.2f]  rad/s  (mech cap = %.1f)\n', ...
    BOX.wr_lo, BOX.wr_hi, LIM.wr_mech_max);
fprintf('  id  in [%7.2f, %7.2f]  A      (fixed reference)\n', ...
    BOX.id_lo, BOX.id_hi);
fprintf('  Vdc in [%7.1f, %7.1f]  V      (sweep step %.1f V)\n', ...
    BOX.Vdc_lo, BOX.Vdc_hi, dDV);

%% ---- 4. Verification: 15 boundary tests in 3 categories --------------
% Per axis (v, ids_ref, Vdc_ref) we now exercise FIVE distinct regimes:
%   (a) inside the envelope (well within all bounds),
%   (b) just inside the upper bound (critical-upper / boundary edge),
%   (c) just inside the lower bound (critical-lower / boundary edge),
%   (d) outside the upper bound (divergent-upper),
%   (e) outside the lower bound (divergent-lower).
% This produces 15 tests grouped into three reportable categories:
%   IN-RANGE  :  3 tests, none of which should violate any constraint
%   CRITICAL  :  6 tests on the envelope edges (one upper + one lower per axis)
%   DIVERGENT :  6 tests outside the envelope (one above + one below per axis)
%
% The earlier 8-test set covered only "above-the-max" divergence on every
% axis except V_dc; this 15-test set additionally covers
%   - v BELOW cut-in (D2),
%   - V_dc ABOVE the capacitor / IGBT rating (D3),
%   - i_ds BELOW the demagnetisation floor (D6),
% and adds three pure in-range scenarios so each reference axis has every
% one of the five regimes (a)-(e) explicitly verified by simulation.
v_hi   = floor(BOX.v_hi);                              % 13 m/s, near v upper edge
v_lo_c = max(LIM.v_min + 0.5, ceil(LIM.v_min*10)/10);  %  4.5 m/s, just above cut-in
Vdc_hi = BOX.Vdc_hi;                                   % 2374.7 V, V_dc upper edge
Vdc_lo = max(round(BOX.Vdc_lo / 10) * 10, ceil(LIM.Vdc_min/10)*10);   % ~ 1010 V floor
ids_lo = LIM.ids_min;                                  % 86 A, demag floor
ids_hi = LIM.ids_max;                                  % 425 A, saturation cap

% Out-of-envelope counter-examples (one above + one below per axis):
v_xs_hi  = 20;     % above electrical envelope 12.45 m/s, inside IEC cut-out 25
v_xs_lo  = 2;      % below cut-in 4 m/s
Vdc_xs_hi = 2700;  % above capacitor rating 2400 V
Vdc_xs_lo = 800;   % below SVPWM floor 980 V
ids_xs_hi = 600;   % above saturation cap 521 A
ids_xs_lo = 50;    % below demag floor 86 A

steady = @(vc) @(t) vc * ones(size(t));

% Note on Test C6 (ids = 86 A) wind speed: at i_ds = 86 A (demagnetisation
% floor) the rotor flux lambda_dr = L_m * i_ds is at its minimum, so the
% torque balance T_e = 1.5*(L_m/L_r)*P*lambda_dr*i_qs requires a much
% larger |i_qs| to deliver the same Cp(v) * P_wind.  Hand calculation:
%   |T_e| = 1.957 * (i_ds/86) * |i_qs|,
%   |T_e| = P_tur(v)/wr_MPPT(v) -> |i_qs| ~ v^2 .
% The |i_s| <= I_s,max constraint pushes the slice-wise feasible v upper
% limit at i_ds = 86 A down to about 7-8 m/s.  We use steady v = 7 m/s
% for C6 to keep the joint (v, i_ds, V_dc) triple inside the 3-D envelope.
TESTS = struct( ...
  'tag', { ...
    'IR1-Nominal',   'IR2-LowMid',    'IR3-HighMid', ...
    'C1-vHi',        'C2-vLo',        'C3-VdcHi',     'C4-VdcLo',     'C5-idsHi',     'C6-idsLo', ...
    'D1-vAbove',     'D2-vBelow',     'D3-VdcAbove',  'D4-VdcBelow',  'D5-idsAbove',  'D6-idsBelow' ...
  }, ...
  'category', { ...
    1, 1, 1, ...
    2, 2, 2, 2, 2, 2, ...
    3, 3, 3, 3, 3, 3 ...
  }, ...
  'desc', { ...
    'paper nominal: v=11 m/s, i_{ds}^{ref}=286 A, V_{dc}^{ref}=1600 V', ...
    'low-mid load: v=8 m/s, i_{ds}^{ref}=200 A, V_{dc}^{ref}=1400 V', ...
    'high-mid load: v=10 m/s, i_{ds}^{ref}=350 A, V_{dc}^{ref}=1900 V', ...
    sprintf('CRIT v upper edge: v=8 -> %d m/s, refs nominal',v_hi), ...
    sprintf('CRIT v lower edge: v=%.1f m/s (just above cut-in %d)',v_lo_c,LIM.v_min), ...
    sprintf('CRIT V_{dc} upper edge: V_{dc}^{ref}=%.0f V (envelope cap)',Vdc_hi), ...
    sprintf('CRIT V_{dc} lower edge: V_{dc}^{ref}=%.0f V (just above SVPWM floor)',Vdc_lo), ...
    sprintf('CRIT i_{ds} upper edge: i_{ds}^{ref}=%.0f A (saturation cap)',ids_hi), ...
    sprintf('CRIT i_{ds} lower edge: i_{ds}^{ref}=%.0f A @ v=7 m/s (demag floor, coupled)',ids_lo), ...
    sprintf('DIV v ABOVE max: v=%d m/s (>%.2f electrical, <%d IEC cut-out)',v_xs_hi,12.45,LIM.v_max), ...
    sprintf('DIV v BELOW min: v=%d m/s (<%d cut-in)',v_xs_lo,LIM.v_min), ...
    sprintf('DIV V_{dc} ABOVE max: V_{dc}^{ref}=%d V (>%d capacitor cap)',Vdc_xs_hi,LIM.Vdc_max), ...
    sprintf('DIV V_{dc} BELOW min: V_{dc}^{ref}=%d V (<%.0f SVPWM floor)',Vdc_xs_lo,LIM.Vdc_min), ...
    sprintf('DIV i_{ds} ABOVE max: i_{ds}^{ref}=%d A (>%d saturation)',ids_xs_hi,LIM.ids_max), ...
    sprintf('DIV i_{ds} BELOW min: i_{ds}^{ref}=%d A (<%d demag floor)',ids_xs_lo,round(LIM.ids_min)) ...
  }, ...
  'ids_ref', { ...
    286, 200, 350, ...
    286, 286, 286, 286, ids_hi, ids_lo, ...
    286, 286, 286, 286, ids_xs_hi, ids_xs_lo ...
  }, ...
  'Vdc_ref', { ...
    1600, 1400, 1900, ...
    1600, 1600, Vdc_hi, Vdc_lo, 1600, 1600, ...
    1600, 1600, Vdc_xs_hi, Vdc_xs_lo, 1600, 1600 ...
  }, ...
  'v_steady', { ...
    11, 8, 10, ...
    v_hi, v_lo_c, 11, 11, 11, 7, ...
    v_xs_hi, v_xs_lo, 11, 11, 11, 11 ...
  }, ...
  'wind', { ...
    steady(11), steady(8), steady(10), ...
    @(t) 8*(t<5)+v_hi*(t>=5), steady(v_lo_c), steady(11), steady(11), steady(11), steady(7), ...
    @(t) 11*(t<5)+v_xs_hi*(t>=5), steady(v_xs_lo), steady(11), steady(11), steady(11), steady(11) ...
  }, ...
  'expect', { ...
    'feasible','feasible','feasible', ...
    'feasible','feasible','feasible','feasible','feasible','feasible', ...
    'INFEASIBLE','INFEASIBLE','INFEASIBLE','INFEASIBLE','INFEASIBLE','INFEASIBLE' ...
  } );

nTests  = numel(TESTS);
T_all   = cell(nTests,1);
X_all   = cell(nTests,1);
verdict = cell(nTests,1);
report  = cell(nTests,1);   % final-state metrics struct

for k = 1:nTests
    p_k = p;
    p_k.ids_ref      = TESTS(k).ids_ref;
    p_k.Vdc_ref      = TESTS(k).Vdc_ref;
    p_k.wind_profile = TESTS(k).wind;
    fprintf('\n=== Verification simulation #%d  [%s] (cat=%d) ===\n', ...
        k, TESTS(k).tag, TESTS(k).category);
    fprintf('  %s   (expect: %s)\n', TESTS(k).desc, TESTS(k).expect);
    [verdict{k}, report{k}, T_all{k}, X_all{k}] = ...
        runAndCheck(p_k, LIM, TESTS(k).v_steady);
    fprintf('  -> verdict: %s\n', verdict{k});
end

%% Summary table
fprintf('\n=== Boundary-point verification summary (15 tests) ===\n');
catLabels = {'IN-RANGE', 'CRITICAL', 'DIVERGENT'};
fprintf('  k tag           cat        ids_f    wr_f    Vdc_f    |i_s|f   |i_g|f   ||m_s||  ||m_g|| verdict\n');
for k = 1:nTests
    r = report{k};
    if isempty(r.t)
        fprintf('  %2d %-13s %-9s   --       --      --       --       --       --       --      %s\n', ...
            k, TESTS(k).tag, catLabels{TESTS(k).category}, verdict{k});
    else
        fprintf('  %2d %-13s %-9s %7.1f %7.2f %8.1f %8.1f %8.1f %8.3f %8.3f  %s\n', ...
            k, TESTS(k).tag, catLabels{TESTS(k).category}, r.ids_f, r.wr_f, r.Vdc_f, ...
            r.is_f, r.ig_f, r.ms_f, r.mg_f, verdict{k});
    end
end

%% ---- 5. Plots (zh + en) ----------------------------------------------
plotDir = pwd;

tit_zh = struct( ...
    'feasmap',  '可行操作包络（概念 A，物理 + 控制器约束求交）', ...
    'margins',  '球面约束与电流裕度（数值越大越宽松）', ...
    'verify_inrange',  '范围内场景（3 张子图）：完全位于可行域内部', ...
    'verify_critical', '临界值场景（6 张子图）：每条参考轴的上下边界各一', ...
    'verify_divergent','发散场景（6 张子图）：每条参考轴的超上限/低于下限各一');
tit_en = struct( ...
    'feasmap',  'Feasibility operating envelope (Concept A)', ...
    'margins',  'Sphere & current margins (larger = more headroom)', ...
    'verify_inrange',  'In-range scenarios (3 panels): all references well inside the envelope', ...
    'verify_critical', 'Critical-edge scenarios (6 panels): upper and lower edge per reference axis', ...
    'verify_divergent','Divergent scenarios (6 panels): above-max and below-min per reference axis');

ax_zh = struct('v','风速 v [m/s]', 'V','直流参考 V_{dc,ref} [V]', ...
               't','时间 t [s]', ...
               'lg_phys','物理范围', 'lg_sphere','球面约束', ...
               'lg_curr','电流额定', 'lg_feas','可行', ...
               'lg_v','v 极限', 'lg_Vdc','V_{dc} 极限');
ax_en = struct('v','wind speed v [m/s]', 'V','DC reference V_{dc,ref} [V]', ...
               't','time t [s]', ...
               'lg_phys','physical box', 'lg_sphere','sphere bind', ...
               'lg_curr','current bind', 'lg_feas','feasible', ...
               'lg_v','v limits', 'lg_Vdc','V_{dc} limits');

% --- Figure 6 : feasibility map ---
saveFeasMap(v_grid, Vdc_grid, fail_code, LIM, BOX, ...
    tit_zh.feasmap, ax_zh, fullfile(plotDir,'fig_feas_envelope_zh.png'));
saveFeasMap(v_grid, Vdc_grid, fail_code, LIM, BOX, ...
    tit_en.feasmap, ax_en, fullfile(plotDir,'fig_feas_envelope_en.png'));

% --- Figure 7 : margins (3 subplots: stator sphere, grid sphere, |is|) ---
saveMarginPanel(v_grid, Vdc_grid, ms2_grid, mg2_grid, is_grid, LIM, ...
    tit_zh.margins, ax_zh, fullfile(plotDir,'fig_feas_margins_zh.png'));
saveMarginPanel(v_grid, Vdc_grid, ms2_grid, mg2_grid, is_grid, LIM, ...
    tit_en.margins, ax_en, fullfile(plotDir,'fig_feas_margins_en.png'));

% --- Figures 8a/8b/8c : boundary-point verification by category ---
% Split the verification plots into three categorised figures so each
% figure visualises a coherent regime.
catIdx = arrayfun(@(t) t.category, TESTS);
idx_in   = find(catIdx == 1);
idx_crit = find(catIdx == 2);
idx_div  = find(catIdx == 3);

saveVerifyPlotCat(T_all(idx_in),   X_all(idx_in),   TESTS(idx_in),   verdict(idx_in),   LIM, ...
    [1 3], tit_zh.verify_inrange,  ax_zh, fullfile(plotDir,'fig_feas_verify_inrange_zh.png'));
saveVerifyPlotCat(T_all(idx_in),   X_all(idx_in),   TESTS(idx_in),   verdict(idx_in),   LIM, ...
    [1 3], tit_en.verify_inrange,  ax_en, fullfile(plotDir,'fig_feas_verify_inrange_en.png'));

saveVerifyPlotCat(T_all(idx_crit), X_all(idx_crit), TESTS(idx_crit), verdict(idx_crit), LIM, ...
    [2 3], tit_zh.verify_critical, ax_zh, fullfile(plotDir,'fig_feas_verify_critical_zh.png'));
saveVerifyPlotCat(T_all(idx_crit), X_all(idx_crit), TESTS(idx_crit), verdict(idx_crit), LIM, ...
    [2 3], tit_en.verify_critical, ax_en, fullfile(plotDir,'fig_feas_verify_critical_en.png'));

saveVerifyPlotCat(T_all(idx_div),  X_all(idx_div),  TESTS(idx_div),  verdict(idx_div),  LIM, ...
    [2 3], tit_zh.verify_divergent,ax_zh, fullfile(plotDir,'fig_feas_verify_divergent_zh.png'));
saveVerifyPlotCat(T_all(idx_div),  X_all(idx_div),  TESTS(idx_div),  verdict(idx_div),  LIM, ...
    [2 3], tit_en.verify_divergent,ax_en, fullfile(plotDir,'fig_feas_verify_divergent_en.png'));

fprintf('\nFigures saved (Chinese / English):\n');
for nm = {'feas_envelope', 'feas_margins', ...
          'feas_verify_inrange', 'feas_verify_critical', 'feas_verify_divergent'}
    fprintf('  fig_%s_zh.png   fig_%s_en.png\n', nm{1}, nm{1});
end

% Persist results
save(fullfile(plotDir, 'feasibility_envelope_results.mat'), ...
    'LIM', 'BOX', 'v_grid', 'Vdc_grid', 'feas', 'fail_code', ...
    'ms2_grid', 'mg2_grid', 'is_grid');

%% =====================================================================
%%                          local functions
%% =====================================================================

function [verdict, info, t, X] = runAndCheck(p, LIM, v_steady)
% Run the 14-state nonlinear simulation under the supplied parameter set
% and grade the resulting trajectory against (a) hardware ratings in LIM
% acting on the actual trajectory metrics, and (b) explicit
% reference-axis violations (since the averaged ODE model does not
% enforce reference-level constraints internally).
% The verdict string is one of:
%   PASS                     - all final-window metrics within ratings
%                              AND no reference-axis violation
%   FAIL: ...                - either a hardware rating is exceeded by
%                              the trajectory, or one or more references
%                              violate the physical envelope
%   DIVERGED                 - trajectory contains Inf/NaN
%   ABORT: <error msg>       - solver threw an exception
%
% v_steady is the steady-state wind speed used for the post-transient
% window (we use it to flag wind-axis violations independent of the
% time-varying wind profile during transients).
    info = struct('t',[],'X',[],'ids_f',NaN,'iqs_f',NaN,'wr_f',NaN, ...
                  'iq_f',NaN,'Vdc_f',NaN,'is_f',NaN,'ig_f',NaN, ...
                  'ms_f',NaN,'mg_f',NaN);
    t = []; X = [];
    try
        [t, X] = runNonlinearSim(p);
    catch ME
        verdict = ['ABORT: ' ME.message];
        return
    end
    if any(~isfinite(X(:)))
        verdict = 'DIVERGED';
        info.t = t; info.X = X;
        return
    end
    n = numel(t); win = max(1, round(n*0.7)):n;
    info.t      = t; info.X = X;
    info.ids_f  = mean(X(win, 1));
    info.iqs_f  = mean(X(win, 2));
    info.wr_f   = mean(X(win, 5));
    info.iq_f   = mean(X(win, 7));
    info.Vdc_f  = mean(X(win, 8));
    info.is_f   = sqrt(info.ids_f^2 + info.iqs_f^2);
    info.ig_f   = abs(info.iq_f);
    info.ms_f   = sqrt(mean(X(win, 9))^2 + mean(X(win,10))^2);
    info.mg_f   = sqrt(mean(X(win,12))^2 + mean(X(win,13))^2);

    flags = {};
    % --- (a) hardware-rating violations on the actual trajectory ---
    if info.is_f  > LIM.is_max
        flags{end+1} = sprintf('|i_s|=%.0f>%d', info.is_f, LIM.is_max);
    end
    if info.ig_f  > LIM.ig_max
        flags{end+1} = sprintf('|i_g|=%.0f>%d', info.ig_f, LIM.ig_max);
    end
    if info.ms_f  > 1
        flags{end+1} = sprintf('||m_s||=%.2f>1', info.ms_f);
    end
    if info.mg_f  > 1
        flags{end+1} = sprintf('||m_g||=%.2f>1', info.mg_f);
    end
    % --- (b) reference-axis (envelope-edge) violations ---
    %   The averaged ODE has no shutdown / saturation / SVPWM clipping,
    %   so we must explicitly flag references that fall outside the
    %   physical envelope of Table 2.
    if v_steady < LIM.v_min
        flags{end+1} = sprintf('v=%.1f<%d (cut-in)', v_steady, LIM.v_min);
    end
    if v_steady > LIM.v_max
        flags{end+1} = sprintf('v=%.1f>%d (cut-out)', v_steady, LIM.v_max);
    end
    if p.Vdc_ref < LIM.Vdc_min
        flags{end+1} = sprintf('V_{dc}^{ref}=%.0f<%.0f (SVPWM floor)', ...
            p.Vdc_ref, LIM.Vdc_min);
    end
    if p.Vdc_ref > LIM.Vdc_max
        flags{end+1} = sprintf('V_{dc}^{ref}=%.0f>%.0f (cap rating)', ...
            p.Vdc_ref, LIM.Vdc_max);
    end
    if p.ids_ref < LIM.ids_min
        flags{end+1} = sprintf('i_{ds}^{ref}=%.0f<%.0f (demag floor)', ...
            p.ids_ref, LIM.ids_min);
    end
    if p.ids_ref > LIM.ids_max
        flags{end+1} = sprintf('i_{ds}^{ref}=%.0f>%d (saturation cap)', ...
            p.ids_ref, LIM.ids_max);
    end
    if isempty(flags)
        verdict = 'PASS';
    else
        verdict = ['FAIL: ' strjoin(flags, '; ')];
    end
end

function [eqv, ms2, mg2, ok_real] = computeEquilibriumDetail(v, p)
% Same as computeEquilibrium in bounds_analysis.m, plus:
%   ms2 = mds^2 + mqs^2  (stator sphere norm squared)
%   mg2 = mdg^2 + mqg^2  (grid sphere norm squared)
%   ok_real = false if quadratic discriminant < 0 (eq. infeasible)
    ok_real = true;
    ms2 = NaN; mg2 = NaN; eqv = nan(14,1);

    wr   = p.lambda_opt * v * p.ng / p.Rb;
    ldr  = p.Lm * p.ids_ref;
    lqr  = 0;
    Vdc  = p.Vdc_ref;
    ids  = p.ids_ref;
    id   = 0;

    lam_i_opt = 1 / (p.lambda_opt + 0.08*p.beta) - 0.035 / (p.beta^3 + 1);
    Cp_opt    = 0.22 * (116 * lam_i_opt - 0.4*p.beta - 5) * exp(-12.5 * lam_i_opt);
    Pwind     = 0.5 * p.rho * pi * p.Rb^2 * v^3;
    Ptur      = Cp_opt * Pwind;
    if abs(wr) < 1e-6, ok_real = false; return, end
    Tm        = -Ptur / wr;
    Te        = p.b * wr + Tm;
    iqs       = Te / (1.5 * (p.Lm / p.Lr) * p.P * ldr);

    aRs  = p.Rr * p.Lm^2 / p.Lr^2 + p.Rs;
    wsl  = (1 / (p.tau_r * p.ids_ref)) * iqs;
    we   = p.P * wr + wsl;
    mds  = (aRs*ids - p.sig*we*iqs - (p.Rr*p.Lm/p.Lr^2)*ldr) / (2*Vdc);
    mqs  = (aRs*iqs + p.sig*we*ids + (p.Lm/p.Lr)*p.P*ldr*wr) / (2*Vdc);
    ms2  = mds^2 + mqs^2;

    S_gen = mds * ids + mqs * iqs;
    aq = p.Rg;
    bq = -p.Um;
    cq = 2*Vdc*S_gen + 2*Vdc^2 / (3*p.Rdc);
    disc = bq^2 - 4*aq*cq;
    if disc < 0
        ok_real = false; return
    end
    iq  = (-bq - sqrt(disc)) / (2*aq);
    mqg = (p.Um - p.Rg*iq) / (2*Vdc);
    mdg = (p.ws*p.Lg*iq) / (2*Vdc);
    mg2 = mdg^2 + mqg^2;
    if ms2 < 0 || mg2 < 0
        ok_real = false; return
    end
    z3  = sqrt(max(0, 1 - ms2));
    z6  = sqrt(max(0, 1 - mg2));

    eqv = [ids; iqs; ldr; lqr; wr; id; iq; Vdc; ...
           mds; mqs; z3; mdg; mqg; z6];
end

function saveFeasMap(vG, VdcG, fc, LIM, ~, ttl, ax, fname)
    fig = figure('Position',[80 80 900 700],'Visible','off','Color','w');
    imagesc(vG, VdcG, fc); axis xy; hold on;
    cmap = [1 1 1;          % 0 feasible
            0.95 0.95 0.95; % 1 wind range
            0.92 0.92 0.92; % 2 Vdc range
            1.0  0.7 0.7;   % 3 stator sphere
            1.0  0.85 0.6;  % 4 grid sphere
            0.7  0.85 1.0;  % 5 |is|
            0.7  1.0  0.7;  % 6 |ig|
            0.85 0.6  1.0]; % 7 rotor over-speed
    colormap(cmap);
    clim([0 7]);

    % Outline feasibility region (where fc==0)
    feasMask = (fc == 0);
    contour(vG, VdcG, double(feasMask), [0.5 0.5], 'k', 'LineWidth', 2.0);

    % Hardware boxes
    plot([LIM.v_min LIM.v_min], ylim, 'b--', 'LineWidth', 1.0);
    plot([LIM.v_max LIM.v_max], ylim, 'b--', 'LineWidth', 1.0);
    plot(xlim, [LIM.Vdc_min LIM.Vdc_min], 'g--', 'LineWidth', 1.0);
    plot(xlim, [LIM.Vdc_max LIM.Vdc_max], 'g--', 'LineWidth', 1.0);

    plot(11, 1600, 'kp', 'MarkerSize', 14, 'MarkerFaceColor','y');
    text(11.2, 1620, 'nominal', 'FontSize', 9);

    xlabel(ax.v); ylabel(ax.V);
    title(ttl);
    grid on;
    cb = colorbar('Ticks', 0:7, 'TickLabels', ...
        {'feasible','wind','Vdc','spr-s','spr-g','|is|','|ig|','wr_{mech}'});
    cb.Label.String = 'binding constraint';

    exportgraphics(fig, fname, 'Resolution', 300);
    close(fig);
end

function saveMarginPanel(vG, VdcG, ms2, mg2, is, LIM, ttl, ax, fname)
    fig = figure('Position',[80 80 1500 450],'Visible','off','Color','w');

    subplot(1,3,1);
    contourf(vG, VdcG, ms2, 0:0.1:1.2, 'LineStyle','none');
    colorbar; hold on;
    contour(vG, VdcG, ms2, [1 1], 'r', 'LineWidth', 1.5);
    xlabel(ax.v); ylabel(ax.V);
    title('||m_s||^2 = m_{ds}^2 + m_{qs}^2');
    grid on; clim([0 1.2]);

    subplot(1,3,2);
    contourf(vG, VdcG, mg2, 0:0.1:1.2, 'LineStyle','none');
    colorbar; hold on;
    contour(vG, VdcG, mg2, [1 1], 'r', 'LineWidth', 1.5);
    xlabel(ax.v); ylabel(ax.V);
    title('||m_g||^2 = m_{dg}^2 + m_{qg}^2');
    grid on; clim([0 1.2]);

    subplot(1,3,3);
    contourf(vG, VdcG, is, 'LineStyle','none');
    colorbar; hold on;
    contour(vG, VdcG, is, [LIM.is_max LIM.is_max], 'r', 'LineWidth', 1.5);
    xlabel(ax.v); ylabel(ax.V);
    title(sprintf('|i_s|  (rating %.0f A)', LIM.is_max));
    grid on;

    sgtitle(ttl);
    exportgraphics(fig, fname, 'Resolution', 300);
    close(fig);
end

function saveVerifyPlotCat(T, X, TESTS, verdict, LIM, layout, ttl, ax, fname)
% Categorised verification plot: arrange n tests into a layout(1) x layout(2)
% grid (rows x cols).  Each subplot has two y-axes:
%   left  axis (black, solid):  |i_s|(t)  with red dashed rating line
%   right axis (blue, dashed):  V_{dc}(t)
% Title shows the test description and the verdict (PASS / FAIL ... /
% DIVERGED / ABORT ...).  Used by Section 9 of bounds_theory.tex to
% generate three figures, one per test category (IN-RANGE / CRITICAL /
% DIVERGENT).
%
% NOTE on axis scaling:  many tests use a steady wind speed, so V_{dc}(t)
% sits essentially flat around its reference value; MATLAB's auto-zoom
% then collapses the right y-axis to an epsilon-thin range and the tick
% labels expose float noise (e.g. 1600.0000000001).  We therefore
% set a per-panel y-range that snaps to the nearest 100 V/100 A and
% enforces a minimum window width, then format ticks as integers.
    nrows = layout(1);
    ncols = layout(2);
    n     = numel(T);
    width = 480 * ncols + 60;
    height = 360 * nrows + 100;
    fig = figure('Position',[40 40 width height],'Visible','off','Color','w');

    for k = 1:n
        ax_h = subplot(nrows, ncols, k); hold on;
        if isempty(X{k})
            text(0.5,0.5,verdict{k},'Units','normalized', ...
                 'HorizontalAlignment','center','Color','r');
            title(sprintf('Test %d  [%s]', k, TESTS(k).tag));
            continue
        end
        is_t = sqrt(X{k}(:,1).^2 + X{k}(:,2).^2);

        % --- Left axis: |i_s|(t) ---
        % Snap to nearest 200 A; minimum span 800 A; always include 0
        % and the rating line LIM.is_max.
        is_top = max([max(is_t)*1.10, LIM.is_max*1.05, 800]);
        is_top = ceil(is_top / 200) * 200;
        yyaxis(ax_h,'left');
        plot(T{k}, is_t, 'k-', 'LineWidth', 1.1);
        yline(LIM.is_max, 'r--', 'LineWidth', 0.9);
        ylabel('|i_s| [A]','Color','k');
        ax_h.YAxis(1).Color = 'k';
        ylim([0, is_top]);
        ax_h.YAxis(1).Exponent = 0;
        ax_h.YAxis(1).TickLabelFormat = '%.0f';

        % --- Right axis: V_{dc}(t) ---
        % Snap to nearest 100 V; minimum span 400 V around mean.
        Vdc_t  = X{k}(:,8);
        v_lo = min(Vdc_t); v_hi = max(Vdc_t); v_mid = 0.5*(v_lo+v_hi);
        v_span_half = max([0.6*(v_hi - v_lo), 200]);
        v_lo_p = floor((v_mid - v_span_half) / 100) * 100;
        v_hi_p = ceil( (v_mid + v_span_half) / 100) * 100;
        yyaxis(ax_h,'right');
        plot(T{k}, Vdc_t, 'b-', 'LineWidth', 1.0);
        ylabel('V_{dc} [V]','Color','b');
        ax_h.YAxis(2).Color = 'b';
        ylim([v_lo_p, v_hi_p]);
        ax_h.YAxis(2).Exponent = 0;
        ax_h.YAxis(2).TickLabelFormat = '%.0f';

        xlabel(ax.t); grid on;
        % Trim long verdict strings so they fit on a single title line.
        v = verdict{k};
        if length(v) > 70, v = [v(1:67) '...']; end
        % Use the test tag as a compact identifier; the description is
        % typically too long to fit in a panel title at a 2x3 layout, so
        % we keep tag + verdict only.
        title({sprintf('%s : %s', TESTS(k).tag, TESTS(k).desc), v}, ...
              'FontSize', 8.5, 'Interpreter','tex');
    end
    sgtitle(ttl, 'FontWeight', 'bold');
    exportgraphics(fig, fname, 'Resolution', 300);
    close(fig);
end

%% --- params and helpers (mirrored from bounds_analysis.m) -------------

function p = setupParams()
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
end

function [t, X] = runNonlinearSim(p)
% Build a 14-state IC at the initial wind speed and integrate
    t_end = 30;
    v0 = p.wind_profile(0);
    p_init = p;
    [eq0, ~, ~, ~] = computeEquilibriumDetail(v0, p_init);
    X0 = eq0;
    odefun = @(t, x) scigRhsXVT(x, p.wind_profile(t), p);
    opts = odeset('RelTol',1e-5,'AbsTol',1e-6,'MaxStep',0.05);
    sol = ode23tb(odefun, [0 t_end], X0, opts);
    t = linspace(0, t_end, 2000)';
    X = deval(sol, t)';
end

function dx = scigRhsXVT(x, v, p)
% Identical to scigRhsXV in bounds_analysis.m
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

    dids = (-aRs*ids + sig*we*iqs + (p.Rr*p.Lm/p.Lr^2)*ldr ...
            + aLm*p.P*lqr*wr + 2*mds*Vdc) / sig;
    diqs = (-aRs*iqs - sig*we*ids + (p.Rr*p.Lm/p.Lr^2)*lqr ...
            - aLm*p.P*ldr*wr + 2*mqs*Vdc) / sig;
    dldr = (p.Rr*p.Lm/p.Lr) * ids - (p.Rr/p.Lr) * ldr ...
           + (we - p.P*wr) * lqr;
    dlqr = (p.Rr*p.Lm/p.Lr) * iqs - (p.Rr/p.Lr) * lqr ...
           - (we - p.P*wr) * ldr;
    Te   = 1.5 * aLm * p.P * (ldr*iqs - lqr*ids);
    dwr  = (Te - p.b*wr - Tm) / p.Jwt;

    Vd = 0;  Vq = p.Um;
    did  = (-p.Rg*id + p.ws*p.Lg*iq - 2*mdg*Vdc + Vd) / p.Lg;
    diq  = (-p.Rg*iq - p.ws*p.Lg*id - 2*mqg*Vdc + Vq) / p.Lg;
    dVdc = (3*(mdg*id + mqg*iq - mds*ids - mqs*iqs) - Vdc/p.Rdc) / p.C;

    e_ids = ids - p.ids_ref;
    e_wr  = wr  - wr_ref;
    s1    = z1^2 + z2^2 + z3^2 - 1;
    dz1   = -p.k1 * e_ids * z3;
    dz2   = -p.k2 * e_wr  * z3;
    dz3   =  p.k1*e_ids*z1 + p.k2*e_wr*z2 - p.c1*s1*z3;

    e_id  = id;
    e_Vdc = Vdc - p.Vdc_ref;
    s2    = z4^2 + z5^2 + z6^2 - 1;
    dz4   = -p.k3 * e_id  * z6;
    dz5   = -p.k4 * e_Vdc * z6;
    dz6   =  p.k3*e_id*z4 + p.k4*e_Vdc*z5 - p.c2*s2*z6;

    dx = [dids; diqs; dldr; dlqr; dwr; did; diq; dVdc; ...
          dz1; dz2; dz3; dz4; dz5; dz6];
end

function Tm = turbineTorque(v, wr, p)
    v       = max(v, 0.1);
    wr_ls   = wr / p.ng;
    lambda  = wr_ls * p.Rb / v;
    lam_i   = 1 / (lambda + 0.08*p.beta) - 0.035 / (p.beta^3 + 1);
    Cp      = 0.22 * (116 * lam_i - 0.4*p.beta - 5) * exp(-12.5 * lam_i);
    Cp      = max(Cp, 0);
    A       = pi * p.Rb^2;
    Pwind   = 0.5 * p.rho * A * v^3;
    Ptur    = Cp * Pwind;
    wr_safe = max(abs(wr), 1e-3) * sign(wr + (wr == 0));
    Tm      = -Ptur / wr_safe;
end
