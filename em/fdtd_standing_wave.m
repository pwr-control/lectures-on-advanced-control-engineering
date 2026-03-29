% =========================================================================
%  FDTD 1D  –  Onda piana in dielettrico puro con parete PEC
%  Schema di Yee leapfrog (staggering in spazio e tempo)
%
%  Equazione governante (sigma = 0):
%       d²E/dx² = (1/c²) * d²E/dt²
%
%  Dominio:   0 <= x <= L
%  BC sinistra:  soft source  E(0,t) = E0 * sin(2*pi*f * t)
%  BC destra:    PEC  -->  E(L,t) = 0   (parete metallica)
%
%  Schema leapfrog (schema di Yee 1D, mezzo non magnetico mu_r=1):
%       H^{n+1/2}_{i+1/2} = H^{n-1/2}_{i+1/2} - (dt/mu0/dx)*(E^n_{i+1} - E^n_i)
%       E^{n+1}_i         = E^n_i               - (dt/eps/dx)*(H^{n+1/2}_{i+1/2} - H^{n+1/2}_{i-1/2})
%
%  Stabilità (CFL):  c * dt / dx <= 1   (qui scelto = 0.99)
%
%  Autore:  generato per Davide Bagnara
%  Data:    2026
% =========================================================================

clear; close all; clc;

%% ── 1. PARAMETRI FISICI ─────────────────────────────────────────────────

c0      = 3e8;                      % velocità della luce [m/s]
mu0     = 4*pi*1e-7;                % permeabilità del vuoto [H/m]
eps0    = 1/(mu0*c0^2);             % permittività del vuoto [F/m]

eps_r   = 4.0;                      % permittività relativa del mezzo (es. PTFE ~ 2.1, vetro ~ 4)
mu_r    = 1.0;                      % permeabilità relativa (dielettrico non magnetico)

eps     = eps_r * eps0;
mu      = mu_r * mu0;
c       = c0 / sqrt(eps_r * mu_r);   % velocità di fase nel mezzo

%% ── 2. GEOMETRIA E SORGENTE ─────────────────────────────────────────────

L       = 1.0;          % lunghezza del dominio [m]
f       = 300e6;        % frequenza della sorgente [Hz]
lambda  = c / f;        % lunghezza d'onda nel mezzo [m]
E0      = 1.0;          % ampiezza del campo E [V/m]

fprintf('Frequenza      : %.1f MHz\n', f/1e6);
fprintf('Lambda nel mezzo: %.4f m\n', lambda);
fprintf('L / lambda     : %.2f\n', L/lambda);
fprintf('Nodi stazionari attesi a x/L =');
n_nodes = floor(2*L/lambda);  % numero di nodi (incluso PEC)
for k = 1:n_nodes
    fprintf(' %.3f', k*lambda/(2));
end
fprintf('\n\n');

%% ── 3. GRIGLIA FDTD ─────────────────────────────────────────────────────
%
%  Griglia di Yee 1D:
%
%   E : nodi interi  -->  x_E(i) = (i-1)*dx,   i = 1 .. Nx+1
%   H : nodi mezzi   -->  x_H(j) = (j-0.5)*dx, j = 1 .. Nx
%
%   |--E--|--H--|--E--|--H--|--E--|  ...  |--E--|--H--|--E--|
%   i=1  j=1  i=2  j=2  i=3             i=Nx  j=Nx i=Nx+1
%
%   BC:  E(1) = sorgente soft source
%        E(Nx+1) = 0  (PEC)
%

CFL     = 0.99;                     % numero di Courant (< 1 per stabilità)
Nx      = 200;                      % celle spaziali
dx      = L / Nx;                   % passo spaziale [m]
dt      = CFL * dx / c;             % passo temporale [s]

% Coordinate per plot
x_E = (0:Nx)  * dx;                 % posizioni nodi E  [m]
x_H = (0.5:1:Nx-0.5) * dx;         % posizioni nodi H  [m]

%% ── 4. COEFFICIENTI AGGIORNAMENTO ───────────────────────────────────────

Ce = dt / (eps * dx);               % coeff. aggiornamento E
Ch = dt / (mu  * dx);               % coeff. aggiornamento H

%% ── 5. INIZIALIZZAZIONE CAMPI ───────────────────────────────────────────

E  = zeros(Nx+1, 1);                % E_{x} nei nodi interi
H  = zeros(Nx,   1);                % H_{y} nei nodi mezzi  (H^{-1/2} = 0)

%% ── 6. TEMPO DI SIMULAZIONE ─────────────────────────────────────────────

T_sim   = 12 / f;                   % simulo 12 periodi (regime stazionario)
Nt      = ceil(T_sim / dt);
t_vec   = (0:Nt-1) * dt;

%% ── 7. STORAGE PER ANIMAZIONE E ANALISI ────────────────────────────────

save_every  = max(1, floor(Nt/300));        % salva ~300 frame
E_history   = zeros(Nx+1, floor(Nt/save_every)+1);
t_history   = zeros(1,    floor(Nt/save_every)+1);
frame_idx   = 1;

% Monitor in un punto fisso (x = L/4)
x_monitor   = L/4;
i_monitor   = round(x_monitor/dx) + 1;
E_monitor   = zeros(1, Nt);

%% ── 8. LOOP FDTD ────────────────────────────────────────────────────────

fprintf('Avvio simulazione FDTD: %d passi temporali...\n', Nt);

for n = 1:Nt

    t_n = (n-1) * dt;

    % ── 8a. Aggiornamento H (passo n-1/2 → n+1/2) ────────────────────
    %   H^{n+1/2}_{j} = H^{n-1/2}_{j} - Ch*(E^n_{j+1} - E^n_{j})
    H = H - Ch * (E(2:end) - E(1:end-1));

    % ── 8b. Aggiornamento E interno (passo n → n+1) ───────────────────
    %   E^{n+1}_{i} = E^n_{i} - Ce*(H^{n+1/2}_{i} - H^{n+1/2}_{i-1})
    E(2:Nx) = E(2:Nx) - Ce * (H(2:end) - H(1:end-1));

    % ── 8c. BC sinistra: soft source ─────────────────────────────────
    %   Aggiunge la sorgente al nodo i=1 (additiva = soft source)
    t_src = t_n + dt;           % tempo dopo aggiornamento E
    E(1) = E(1) + E0 * sin(2*pi*f * t_src);

    % ── 8d. BC destra: PEC ───────────────────────────────────────────
    E(Nx+1) = 0;

    % ── 8e. Monitor ──────────────────────────────────────────────────
    E_monitor(n) = E(i_monitor);

    % ── 8f. Salvataggio frame ─────────────────────────────────────────
    if mod(n, save_every) == 0
        E_history(:, frame_idx) = E;
        t_history(frame_idx)    = t_src;
        frame_idx = frame_idx + 1;
    end
end

fprintf('Simulazione completata.\n\n');

%% ── 9. SOLUZIONE ANALITICA (regime stazionario) ─────────────────────────
%
%  Con PEC in x=L e sorgente in x=0 la soluzione stazionaria è:
%
%    E(x,t) = A * sin(k*(L-x)) * sin(2*pi*f*t)
%
%  dove k = 2*pi/lambda = 2*pi*f/c
%  L'ampiezza A dipende dallo schema soft source (non è banalmente E0),
%  la calcoliamo dal campo simulato nell'ultimo periodo.
%

k_wave = 2*pi*f / c;

% Ampiezza dell'onda stazionaria stimata dal massimo temporale simulato
E_max_num = max(abs(E_history), [], 2);   % profilo |E|_max spaziale

% Ampiezza analitica normalizzata: sin(k*(L-x))
E_analytical_shape = abs(sin(k_wave * (L - x_E')));
% Scala per minimizzare errore quadratico (least squares)
A_fit = (E_analytical_shape' * E_max_num) / (E_analytical_shape' * E_analytical_shape);

fprintf('Ampiezza onda stazionaria (fit): A = %.4f V/m\n', A_fit);
fprintf('k*lambda = %.4f (dovrebbe essere 2*pi = %.4f)\n', k_wave*lambda, 2*pi);

%% ── 10. FIGURA 1: Evoluzione temporale ──────────────────────────────────

fig1 = figure('Name','FDTD – Onda Stazionaria 1D','Color','white','Position',[50 50 1100 700]);

% ── 10a. Animazione campo E ───────────────────────────────────────────────
subplot(3,2,[1 2]);
h_E_line  = plot(x_E, E_history(:,1), 'b-', 'LineWidth', 1.8); hold on;
h_E_anal  = plot(x_E, A_fit*sin(k_wave*(L-x_E))*sin(2*pi*f*t_history(1)), ...
                 'r--', 'LineWidth', 1.4);
h_E_env_p = plot(x_E,  A_fit*abs(sin(k_wave*(L-x_E))), 'k:', 'LineWidth', 1.2);
h_E_env_m = plot(x_E, -A_fit*abs(sin(k_wave*(L-x_E))), 'k:', 'LineWidth', 1.2);
xline(L, 'k-', 'LineWidth', 2.5);
xline(0, 'm--', 'LineWidth', 1.5);
ylabel('E_y  [V/m]');
xlabel('x  [m]');
title('Campo Elettrico E_y(x,t)  –  Onda stazionaria in dielettrico puro','FontSize',11);
legend('FDTD','Analitico','Inviluppo ±','','PEC','Sorgente','Location','northwest');
ylim([-2.5*A_fit, 2.5*A_fit]);
xlim([0 L]);
grid on; box on;
h_time_txt = text(0.01*L, 2.2*A_fit, '', 'FontSize', 10, 'Color', 'k');

% ── 10b. Monitor temporale ────────────────────────────────────────────────
subplot(3,2,3);
plot(t_vec*1e9, E_monitor, 'b-', 'LineWidth', 1.2);
xlabel('t  [ns]');
ylabel('E_y  [V/m]');
title(sprintf('Monitor E_y in x = %.3f m  (x/\\lambda = %.2f)', x_monitor, x_monitor/lambda));
grid on; box on;

% ── 10c. Profilo |E|_max spaziale vs analitico ───────────────────────────
subplot(3,2,4);
plot(x_E, E_max_num, 'b-', 'LineWidth', 1.8); hold on;
plot(x_E, A_fit * E_analytical_shape, 'r--', 'LineWidth', 1.4);
xlabel('x  [m]');
ylabel('|E|_{max}  [V/m]');
title('Inviluppo spaziale  |E(x)|_{max}');
legend('FDTD','A·|sin(k(L-x))|','Location','northwest');
xline(L,'k-','LineWidth',2); xline(0,'m--','LineWidth',1.5);
grid on; box on; xlim([0 L]);

% ── 10d. Spettro del monitor ──────────────────────────────────────────────
subplot(3,2,5);
Nfft    = 2^nextpow2(Nt);
E_fft   = abs(fft(E_monitor, Nfft));
f_axis  = (0:Nfft/2-1) / (Nfft*dt) / 1e6;   % MHz
E_fft_half = E_fft(1:Nfft/2);
plot(f_axis, 20*log10(E_fft_half/max(E_fft_half)), 'b-', 'LineWidth', 1.2);
xline(f/1e6, 'r--', 'LineWidth', 1.5);
xlabel('Frequenza  [MHz]');
ylabel('|FFT|  [dB]');
title('Spettro FFT del monitor');
xlim([0 3*f/1e6]); ylim([-60 5]);
grid on; box on;

% ── 10e. Errore FDTD vs analitico (regime stazionario) ───────────────────
subplot(3,2,6);
% Usa gli ultimi frame (regime stazionario raggiunto)
n_last   = frame_idx - 1;
n_period = round(1/(f*save_every*dt));
idx_ss   = max(1, n_last - n_period) : n_last;    % ultimi ~ 1 periodo
E_rms_num  = sqrt(mean(E_history(:, idx_ss).^2, 2));
E_rms_anal = A_fit / sqrt(2) * E_analytical_shape;
rel_err    = abs(E_rms_num - E_rms_anal) ./ (A_fit/sqrt(2) + eps);
plot(x_E, rel_err * 100, 'g-', 'LineWidth', 1.5);
xlabel('x  [m]');
ylabel('Errore relativo  [%]');
title(sprintf('Errore FDTD vs analitico  (CFL = %.2f, Nx = %d)', CFL, Nx));
grid on; box on; xlim([0 L]);

%% ── 11. ANIMAZIONE INTERATTIVA ──────────────────────────────────────────

fprintf('Animazione in corso... (chiudi la finestra per fermarla)\n');

n_frames = frame_idx - 1;
for fr = 1:n_frames
    if ~ishandle(fig1), break; end

    t_fr = t_history(fr);

    % Aggiorna linee
    set(h_E_line, 'YData', E_history(:, fr));
    set(h_E_anal, 'YData', A_fit * sin(k_wave*(L-x_E)) * sin(2*pi*f*t_fr));
    set(h_time_txt, 'String', sprintf('t = %.2f ns  (%.1f T)', t_fr*1e9, t_fr*f));

    drawnow limitrate;
    pause(0.005);
end

%% ── 12. FIGURA 2: Mappa spazio-tempo (kymograph) ───────────────────────

figure('Name','Kymograph E(x,t)','Color','white','Position',[200 200 900 450]);

% Seleziona ultimi 2 periodi per il kymograph
n_2T     = round(2/(f*save_every*dt));
idx_kymo = max(1, n_frames-n_2T) : n_frames;

imagesc(x_E, t_history(idx_kymo)*1e9, E_history(:, idx_kymo)');
set(gca,'YDir','normal');
colormap(redblue_colormap());
colorbar;
clim([-A_fit A_fit]);
xlabel('x  [m]','FontSize',12);
ylabel('t  [ns]','FontSize',12);
title('Kymograph  E_y(x,t)  –  ultimi 2 periodi (regime stazionario)','FontSize',12);
xline(L, 'w-', 'LineWidth', 2);
xline(0, 'w--', 'LineWidth', 1.5);

%% ── FUNZIONE LOCALE: colormap divergente rosso-bianco-blu ───────────────

function cmap = redblue_colormap(n)
    if nargin < 1, n = 256; end
    top    = [linspace(1,0,n/2)', linspace(1,0,n/2)', ones(n/2,1)];
    bottom = [ones(n/2,1), linspace(0,1,n/2)', linspace(0,1,n/2)'];
    cmap   = [top; bottom];
end
