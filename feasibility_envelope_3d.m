close all; clear; clc;
p = setupParams();

%% ---- 1. Sweep grid -----------------------------------------------------
% Resolution matched to the 2-D version (feasibility_envelope.m): Δv≈0.35
% m/s, ΔVdc≈30 V.  ids_grid is anchored at the analytical floor (86 A from
% Table 1 / 0.3 * ids_ref) and ceiling (425 A from λ_sat/L_m * 0.82) so
% the scan reports floor/ceiling exactly without snapping artefacts.
v_grid   = linspace(3,   18,  44);    % wind speed [m/s], Δv = 0.349
ids_grid = linspace(86,  425, 35);    % ids_ref    [A], Δi = 9.97, includes 86 & 425 exactly
Vdc_grid = linspace(800, 2500, 57);   % Vdc_ref    [V], ΔVdc = 30.4 (matches 2D)
nV = numel(v_grid); nI = numel(ids_grid); nD = numel(Vdc_grid);

%% ---- 2. Constraint values ---------------------------------------------
% Identical to the 2-D version (feasibility_envelope.m); see the lengthy
% comments there for the engineering source of every value.
LIM = struct();
LIM.v_min   = 4;
LIM.v_max   = 25;
LIM.Vdc_min = sqrt(6) * 400;        % = 979.8 V
LIM.Vdc_max = 2400;
LIM.is_max  = 2920;                 % stator nameplate (PF=0.85, eta=0.95)
LIM.ig_max  = 2357;                 % grid converter (unity PF)
% Absolute ids limits anchored on the *true* machine rated magnetising
% current i_m,rated = V_m/(omega_s*L_m) = 347 A:
%   lambda_sat = 1.5 * V_m/omega_s = 2.702 Wb
%   i_ds,max,sat = lambda_sat / L_m = 521 A
%   adopt 425 = 0.82 * 521 (~18% saturation margin)
% Floor: 0.3 * ids_ref = 86 A (anti-demagnetisation).
% Reference: Krause/Wasynczuk/Sudhoff, Analysis of Electric Machinery, Ch.4.
LIM.ids_abs_min = 0.3 * 286;        % = 86 A
LIM.ids_abs_max = 425;              % = 0.82 * lambda_sat / L_m
% Mechanical over-speed cap.  Full converter decouples generator
% electrical from grid frequency: limit comes from generator/gearbox
% nameplate, NOT from ws/P.  168 rad/s ~ 1.07 * 1500 rpm typical
% nameplate, well below IEC 61400-1 trip threshold ~235-267 rad/s.
LIM.wr_mech_max = 168;              % steady-state mech cap [rad/s]

fprintf('=== 3D feasibility envelope sweep ===\n');
fprintf('  v       grid : %d points in [%.1f, %.1f] m/s\n', nV, v_grid(1), v_grid(end));
fprintf('  ids_ref grid : %d points in [%.0f, %.0f] A\n',   nI, ids_grid(1), ids_grid(end));
fprintf('  Vdc_ref grid : %d points in [%.0f, %.0f] V\n',   nD, Vdc_grid(1), Vdc_grid(end));
fprintf('  total cells  : %d\n\n', nV*nI*nD);

%% ---- 3. Sweep ----------------------------------------------------------
feas3 = false(nD, nV, nI);
fail3 = zeros(nD, nV, nI);
% fail codes (consistent with the 2D version):
% 0 = feasible, 1 = wind range, 2 = Vdc range,
% 3 = stator sphere, 4 = grid sphere,
% 5 = |i_s| > rating, 6 = |i_g| > rating,
% 7 = wr_mech (rotor over-speed),
% 8 = ids floor/ceiling, -1 = numerical infeasibility

for ii = 1:nI
    p_loc = p;
    p_loc.ids_ref = ids_grid(ii);
    for iv = 1:nV
        v = v_grid(iv);
        if v < LIM.v_min || v > LIM.v_max
            fail3(:, iv, ii) = 1;  continue
        end
        for id_ = 1:nD
            Vdc = Vdc_grid(id_);
            p_loc.Vdc_ref = Vdc;

            if Vdc < LIM.Vdc_min || Vdc > LIM.Vdc_max
                fail3(id_, iv, ii) = 2; continue
            end

            try
                [eqv, ms2, mg2, ok_real] = computeEquilibriumDetail(v, p_loc);
            catch
                fail3(id_, iv, ii) = -1; continue
            end
            if ~ok_real
                fail3(id_, iv, ii) = -1; continue
            end

            ids_e = eqv(1); iqs_e = eqv(2); wr_e = eqv(5); iq_e = eqv(7);
            is_mag = sqrt(ids_e^2 + iqs_e^2);
            ig_mag = abs(iq_e);

            if ms2 > 1,                        fail3(id_,iv,ii) = 3; continue, end
            if mg2 > 1,                        fail3(id_,iv,ii) = 4; continue, end
            if is_mag > LIM.is_max,            fail3(id_,iv,ii) = 5; continue, end
            if ig_mag > LIM.ig_max,            fail3(id_,iv,ii) = 6; continue, end
            if wr_e   > LIM.wr_mech_max,       fail3(id_,iv,ii) = 7; continue, end
            if ids_e < LIM.ids_abs_min || ids_e > LIM.ids_abs_max
                fail3(id_,iv,ii) = 8; continue
            end

            feas3(id_, iv, ii) = true;
        end
    end
    if mod(ii, 5) == 0
        fprintf('  ids_ref slice %d / %d done\n', ii, nI);
    end
end

%% ---- 4. Project to bounding box ---------------------------------------
[idxD, idxV, idxI] = ind2sub(size(feas3), find(feas3));
v_feas   = v_grid(unique(idxV));
ids_feas = ids_grid(unique(idxI));
Vdc_feas = Vdc_grid(unique(idxD));
wr_feas  = p.lambda_opt * v_feas * p.ng / p.Rb;

n_total  = numel(feas3);
n_feas   = numel(idxV);

fprintf('\n=== 3D feasibility envelope summary ===\n');
fprintf('Feasible cells : %d / %d  (%.2f%%)\n', n_feas, n_total, 100*n_feas/n_total);
fprintf('  v       feasible in [%.2f, %.2f] m/s\n',  min(v_feas),   max(v_feas));
fprintf('  ids_ref feasible in [%.0f, %.0f] A\n',    min(ids_feas), max(ids_feas));
fprintf('  Vdc_ref feasible in [%.0f, %.0f] V\n',    min(Vdc_feas), max(Vdc_feas));
fprintf('  -> wr   feasible in [%.2f, %.2f] rad/s (= MPPT(v))\n', ...
    min(wr_feas), max(wr_feas));

%% ---- 5. State operating box (Concept A, 3-D) --------------------------
% Sweep step sizes (used to round the box edges to grid resolution).
dv  = v_grid(2)   - v_grid(1);
dI  = ids_grid(2) - ids_grid(1);
dDV = Vdc_grid(2) - Vdc_grid(1);

BOX = struct();
BOX.v_lo   = min(v_feas);    BOX.v_hi   = max(v_feas);
BOX.ids_lo = min(ids_feas);  BOX.ids_hi = max(ids_feas);
BOX.wr_lo  = min(wr_feas);   BOX.wr_hi  = max(wr_feas);
BOX.id_lo  = 0;              BOX.id_hi  = 0;
BOX.Vdc_lo = min(Vdc_feas);  BOX.Vdc_hi = max(Vdc_feas);
fprintf('\n=== State operating box (3D Concept A) ===\n');
fprintf('  v   in [%7.2f, %7.2f] m/s   (sweep step %.2f)\n',  BOX.v_lo,   BOX.v_hi,   dv);
fprintf('  ids in [%7.1f, %7.1f] A     (sweep step %.1f)\n',  BOX.ids_lo, BOX.ids_hi, dI);
fprintf('  wr  in [%7.2f, %7.2f] rad/s (mech cap %.1f)\n',    BOX.wr_lo,  BOX.wr_hi, LIM.wr_mech_max);
fprintf('  id  in [%7.1f, %7.1f] A     (fixed reference)\n',  BOX.id_lo,  BOX.id_hi);
fprintf('  Vdc in [%7.1f, %7.1f] V     (sweep step %.1f)\n',  BOX.Vdc_lo, BOX.Vdc_hi, dDV);

%% ---- 6. Plots (zh + en) ----------------------------------------------
plotDir = pwd;

tit_zh = struct( ...
    'slices',  '可行域多切片（每个 i_{ds}^{ref} 一张子图）', ...
    'volume',  '3-D 可行包络与投影（Concept A 扩展）', ...
    'count',   '可行 (v,V_{dc}^{ref}) 格数随 i_{ds}^{ref} 变化');
tit_en = struct( ...
    'slices',  'Feasibility slices (one panel per i_{ds}^{ref})', ...
    'volume',  '3-D feasibility envelope and projections', ...
    'count',   'feasible (v,V_{dc}^{ref}) cells vs i_{ds}^{ref}');

ax_zh = struct('v','风速 v [m/s]', 'V','直流参考 V_{dc,ref} [V]', ...
               'I','d 轴电流参考 i_{ds}^{ref} [A]', 'cnt','可行格数');
ax_en = struct('v','wind speed v [m/s]', 'V','DC reference V_{dc,ref} [V]', ...
               'I','d-axis ref i_{ds}^{ref} [A]', 'cnt','feasible cells');

% Slices to draw, picking 6 ids_ref values that span the new feasible
% range [86, 425] A imposed by the corrected LIM constants.
slice_targets = [86 150 220 286 360 425];
slice_idx = arrayfun(@(t) find(abs(ids_grid - t) == min(abs(ids_grid - t)), 1), ...
                     slice_targets);

% --- Figure 1 : 6-slice (v, Vdc_ref) feasibility maps -------------------
saveSlicePanel(v_grid, Vdc_grid, ids_grid, fail3, slice_idx, ...
    LIM, tit_zh.slices, ax_zh, fullfile(plotDir, 'fig_feas_3d_slices_zh.png'));
saveSlicePanel(v_grid, Vdc_grid, ids_grid, fail3, slice_idx, ...
    LIM, tit_en.slices, ax_en, fullfile(plotDir, 'fig_feas_3d_slices_en.png'));

% --- Figure 2 : 3-D scatter + bar chart projection ----------------------
saveVolumeFig(v_grid, ids_grid, Vdc_grid, feas3, ...
    tit_zh.volume, tit_zh.count, ax_zh, ...
    fullfile(plotDir, 'fig_feas_3d_volume_zh.png'));
saveVolumeFig(v_grid, ids_grid, Vdc_grid, feas3, ...
    tit_en.volume, tit_en.count, ax_en, ...
    fullfile(plotDir, 'fig_feas_3d_volume_en.png'));

fprintf('\nFigures saved (Chinese / English):\n');
fprintf('  fig_feas_3d_slices_zh.png   fig_feas_3d_slices_en.png\n');
fprintf('  fig_feas_3d_volume_zh.png   fig_feas_3d_volume_en.png\n');

% Persist results
save(fullfile(plotDir, 'feasibility_envelope_3d_results.mat'), ...
    'LIM', 'BOX', 'v_grid', 'ids_grid', 'Vdc_grid', 'feas3', 'fail3');

%% =====================================================================
%%                          local functions
%% =====================================================================

function saveSlicePanel(vG, VdcG, idsG, fail3, slice_idx, LIM, ttl, ax, fname)
    fig = figure('Position', [40 40 1500 800], 'Visible','off', 'Color','w');
    cmap = [1 1 1;          % 0 feasible
            0.95 0.95 0.95; % 1 wind range
            0.92 0.92 0.92; % 2 Vdc range
            1.0  0.7 0.7;   % 3 stator sphere
            1.0  0.85 0.6;  % 4 grid sphere
            0.7  0.85 1.0;  % 5 |is|
            0.7  1.0  0.7;  % 6 |ig|
            0.85 0.6 1.0;   % 7 wr_mech
            0.6  0.6 0.6];  % 8 ids floor/ceil
    nslice = numel(slice_idx);
    for k = 1:nslice
        subplot(2, 3, k);
        ii = slice_idx(k);
        slice = squeeze(fail3(:, :, ii));
        imagesc(vG, VdcG, slice); axis xy; hold on;
        colormap(gca, cmap);
        clim([0 8]);
        contour(vG, VdcG, double(slice == 0), [0.5 0.5], 'k', 'LineWidth', 1.5);
        plot([LIM.v_min LIM.v_min], ylim, 'b--', 'LineWidth', 0.8);
        plot([LIM.v_max LIM.v_max], ylim, 'b--', 'LineWidth', 0.8);
        plot(xlim, [LIM.Vdc_min LIM.Vdc_min], 'g--', 'LineWidth', 0.8);
        plot(xlim, [LIM.Vdc_max LIM.Vdc_max], 'g--', 'LineWidth', 0.8);
        xlabel(ax.v); ylabel(ax.V);
        title(sprintf('i_{ds}^{ref} = %.0f A   feas frac = %.0f%%', ...
            idsG(ii), 100*nnz(slice == 0)/numel(slice)));
        grid on;
    end
    sgtitle(ttl);
    exportgraphics(fig, fname, 'Resolution', 300);
    close(fig);
end

function saveVolumeFig(vG, idsG, VdcG, feas3, ttl, cnt_ttl, ax, fname)
    fig = figure('Position', [40 40 1500 600], 'Visible','off', 'Color','w');

    % --- (a) 3-D scatter of feasible cells ---
    subplot(1, 3, 1);
    [I, V, II] = ind2sub(size(feas3), find(feas3));
    scatter3(vG(V), idsG(II), VdcG(I), 4, VdcG(I), 'filled');
    colormap(gca, 'parula');
    xlabel(ax.v); ylabel(ax.I); zlabel(ax.V);
    title('feasible (v, i_{ds}^{ref}, V_{dc}^{ref}) cells');
    grid on; view(40, 25);

    % --- (b) projection: count per ids_ref ---
    subplot(1, 3, 2);
    n_per_ids = squeeze(sum(sum(feas3, 1), 2));
    bar(idsG, n_per_ids, 'FaceColor', [0.3 0.5 0.85]);
    xlabel(ax.I); ylabel(ax.cnt);
    title(cnt_ttl);
    grid on;

    % --- (c) feasible ids_ref interval (highlight) ---
    subplot(1, 3, 3);
    feas_per_ids = n_per_ids > 0;
    if any(feas_per_ids)
        ids_lo = idsG(find(feas_per_ids, 1, 'first'));
        ids_hi = idsG(find(feas_per_ids, 1, 'last'));
        plot(idsG, double(feas_per_ids), 'o-', 'LineWidth', 1.5, ...
            'MarkerFaceColor','b'); hold on;
        yline(0.5, 'r--');
        xline(ids_lo, 'g--', 'LineWidth', 1.0);
        xline(ids_hi, 'g--', 'LineWidth', 1.0);
        title(sprintf('i_{ds}^{ref} feasible in [%.0f, %.0f] A', ids_lo, ids_hi));
    else
        title('no feasible i_{ds}^{ref}');
    end
    xlabel(ax.I); ylabel('feasible (0/1)');
    grid on; ylim([-0.1, 1.2]);

    sgtitle(ttl);
    exportgraphics(fig, fname, 'Resolution', 300);
    close(fig);
end

%% --- params and equilibrium helpers (mirrored from feasibility_envelope.m) ---

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
    p.k1 = +0.03;  p.k2 = +0.003;
    p.k3 = -1;     p.k4 = -0.01;
    p.c1 = 100;    p.c2 = 100;
end

function [eqv, ms2, mg2, ok_real] = computeEquilibriumDetail(v, p)
% Closed-form equilibrium of the 14-state closed loop, plus
% sphere-norm squared values for both sides.
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

    S_gen = mds*ids + mqs*iqs;
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
