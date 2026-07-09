clc
clearvars
close all

%% ---- SIM CONTROL VALUES ---- %%

numpts = 200; % number of axial nodes
numpts_radial = 120; % number of radial annular sections
Tf = 5; % burn time of engine seconds
numpts_time = 30;

%% ---- Bartz Equation Parameters --- %%

P0_bartz     = 5e5;        % Chamber  pressure [Pa]
gamma_bartz  = 1.285;      % Ratio of specific heats
Cp_gas_bartz = 2811.9;     % Specific heat of gas [J/kg.K]
omega_bartz  = 0.6;        % fixed 
rc_curv      = 0.005;      % Throat radius of curvature [m]
Pr           = 0.51;       % Prandtl number at stagnation conditions
mu = 0.00028743; % viscosity [kg/(m*s)] 
mdot_total = 0.1043; % Total Mass flow rate [Kg/s] 
Twg_bartz = 800;   % Hot-side wall temperature, assumed based on papers[K]

%% ---- Structural Parameters --- %%

k_gr = 16.2; % conductivity of Material [W/(M*K)]
rho_gr = 8000; % density of material [kg/m3]
Cp_gr = 500; % Cp of material
Tinit = 300; % initial wall temperature [K] 


%% --- Define Engine geometry (conical, no fillet) --- %%
Rc = 0.0324;   % chamber radius, [m]
Rt = 0.0102;   % throat radius, [m]
Re = 0.0152;   % exit radius, [m]

Lc     = 0.0909; % chamber length, [m]
L_conv = 0.0475; % converging section length, [m]
L_div  = 0.0185; % diverging section length, [m]

wall_thickness_mm = 10; % radial wall thickness, [mm]

%% --- Number of divisions per section (edit these to change resolution) --- %%

n_chamber = 4; % points in flat chamber section
n_conv    = 8; % points in converging section (ends at throat)
n_div     = 9; % points in diverging section (ends at exit)


%% -------------------------- %%

Le = Lc + L_conv + L_div; % total engine length, [m]

x_chamber_end = Lc;
x_throat      = Lc + L_conv;
x_exit        = Le;
x_chamber = linspace(0, x_chamber_end, n_chamber+1); x_chamber = x_chamber(1:end-1);
x_conv    = linspace(x_chamber_end, x_throat, n_conv+1); x_conv = x_conv(1:end-1);
x_div     = linspace(x_throat, x_exit, n_div);

x = [x_chamber, x_conv, x_div]; % chamber -> exit

interpts = linspace(0,Le,numpts);

% --- Diameter_Array built automatically from Rc, Rt, Re (conical, no fillet) ---
r_chamber = Rc*ones(size(x_chamber));
r_conv    = linspace(Rc, Rt, n_conv+1); r_conv = r_conv(1:end-1);
r_div     = linspace(Rt, Re, n_div);

Diameter_Array = (2*[r_chamber, r_conv, r_div])'; % meters, chamber -> exit, column vector


Spline_Dia = pchip(x,Diameter_Array,interpts);

figure(1)
plot(x,Diameter_Array/2,'k-',LineWidth=1.5)
hold on
plot(x,-Diameter_Array/2,'k-',LineWidth=1.5)
plot(interpts,Spline_Dia/2,'b-',LineWidth=1.5)
plot(interpts,-Spline_Dia/2,'b-',LineWidth=1.5)
xlabel('Axial position (m)')
ylabel('Radius (m)')
title('Thrust Chamber Cross Section (m)')
grid on
axis equal
hold off

Metric_Dia = Spline_Dia;
x_metric = interpts;
Le_Metric = Le;

%% --- Gas parameter interpolation --- %%

x_hga = [0, x_chamber_end, x_throat, x_exit]; 
hga = [1882.76,1881.85,1814,1732]; % recovery temperature at the above mentioned stations
Prga = [0.5084,0.5084,0.5074,0.5074]; % Prandtl number at the above mentioned stations
kga = [0.3402,0.3370,0.3036,0.2227]; % conductivity of gases at the above mentioned stations

%% -------------------------- %%

hotgasarray = pchip(x_hga,hga,x_metric);
pr_gas_array = pchip(x_hga,Prga,x_metric);
k_array = pchip(x_hga,kga,x_metric);


%Initialize and solve for HT coefficient using the Bartz correlation
%(same method as transient_thermal_rad_throat.m), evaluated at each
%axial station along x_metric. Result stored in hg_array_Nusselt to
%keep all downstream variable names unchanged.

R_local   = Metric_Dia/2;           % local radius [m], chamber -> exit
A_local   = pi*R_local.^2;          % local cross-sectional area [m^2]
At_bartz  = pi*Rt^2;                % throat area [m^2]

[~, idx_throat_bartz] = min(A_local);   % throat index (minimum area)

% --- Local Mach number from area-Mach relation (subsonic upstream of
%     throat, supersonic downstream), same as reference file ---
M_array = zeros(size(Metric_Dia));

areaMachFun_bartz = @(M, AR) (1./M) .* ( (2/(gamma_bartz+1)) .* ...
    (1 + (gamma_bartz-1)/2 .* M.^2) ).^((gamma_bartz+1)/(2*(gamma_bartz-1))) - AR;

for i = 1:numpts
    AR = A_local(i)/At_bartz;   % local area ratio

    if i == idx_throat_bartz
        M_array(i) = 1.0;
        continue
    end

    if AR < 1
        AR = 1;   % numerical noise near throat -> clamp to sonic
    end

    if i < idx_throat_bartz
        % Subsonic branch
        if AR <= 1
            M_array(i) = 1.0;
        else
            M_array(i) = fzero(@(M) areaMachFun_bartz(M, AR), [1e-4, 1-1e-6]);
        end
    else
        % Supersonic branch
        if AR <= 1
            M_array(i) = 1.0;
        else
            M_array(i) = fzero(@(M) areaMachFun_bartz(M, AR), [1+1e-6, 15]);
        end
    end
end

% --- Characteristic velocity c* (from chamber conditions) ---
c_star_bartz = P0_bartz * At_bartz / mdot_total;   % [m/s]

% --- Sigma correction factor (uses local gas stagnation temp hotgasarray
%     as T0g, and Twg assumed equal to gas-side array value at hand;
%     here Twg taken as the wall-adjacent reference 800 K used in the
%     reference script) ---


sigma_array = zeros(size(Metric_Dia));
for i = 1:numpts
    Mi = M_array(i);
    T0g_local = hotgasarray(i);   % local hot-gas stagnation temp [K]
    term1 = 0.5*(Twg_bartz/T0g_local)*(1 + (gamma_bartz-1)/2*Mi^2) + 0.5;
    term2 = 1 + (gamma_bartz-1)/2*Mi^2;
    sigma_array(i) = 1 / ( term1^(0.8 - 0.2*omega_bartz) * term2^(0.2*omega_bartz) );
end

% --- Bartz convective heat transfer coefficient ---
D_star_bartz = 2*Rt;   % throat (hydraulic) diameter [m]

const_term_bartz = (0.026 / D_star_bartz^0.2) * (mu^0.2 * Cp_gas_bartz / Pr^0.6) * ...
    (P0_bartz / c_star_bartz)^0.8 * (D_star_bartz / rc_curv)^0.1;

hg_array_Nusselt = zeros(size(Metric_Dia));
for i = 1:numpts
    Prantl = pr_gas_array(i);   % retained for compatibility (unused in Bartz eq itself)
    k = k_array(i);              % retained for compatibility (unused in Bartz eq itself)
    area_ratio_term_bartz = (At_bartz / A_local(i))^0.9;
    hg_array_Nusselt(i) = const_term_bartz * area_ratio_term_bartz * sigma_array(i);
end
figure(3)
plot(x_metric,k_array,'LineWidth',1.5)

grid on

figure(2)
plot(x_metric,hg_array_Nusselt/1000,'r-',Linewidth=2) 
xlabel('Axial position (m)','FontSize',16)
ylabel('Convective Coefficient (kW/(m^2K))','FontSize',16)        
grid on
title('Bartz Coefficient vs Axial Position','FontSize',16)
hold on
[hg_peak, idx_hg_peak] = max(hg_array_Nusselt);
plot(x_metric(idx_hg_peak), hg_peak/1000, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
text(x_metric(idx_hg_peak), hg_peak/1000, ...
    sprintf('  Peak = %.2f kW/(m^2K) @ x = %.4f m', hg_peak/1000, x_metric(idx_hg_peak)), ...
    'FontSize',12,'VerticalAlignment','bottom')
hold off

% Outer wall now traces the inner contour at a constant radial thickness
% (curved jacket), instead of a flat cylinder.
wall_thickness_m = wall_thickness_mm/1000; % m
Outer_Dia = Metric_Dia + 2*wall_thickness_m; % Metric_Dia already in meters

Tg_array = ones(size(x_metric))*400;



timepoints = linspace(0,Tf,numpts_time);

Tdiskarray = zeros([numpts,numpts_time]);


tau_array = ones(size(x_metric));
biot_array = ones(size(x_metric));
for i = 1:numpts
    hotgas_temp = hotgasarray(i);
    hg = hg_array_Nusselt(i);
    india = Metric_Dia(i);
    pakistan = Outer_Dia(i);
    V = (Le_Metric/numpts)*(pi/4)*(pakistan^2 - india^2);
    tau = (hg*(pi*india*(Le_Metric/numpts)))/(rho_gr*Cp_gr*V);
    biot_array(i) = (hg*V)/(k_gr*(pi*india*(Le_Metric/numpts)));
    tau_array(i) = tau;
    for j = 1:numpts_time
        Tdiskarray(i,j) = hotgas_temp + (297 - hotgas_temp)*exp(-tau*timepoints(j));

    end
end


%%%%%%%%%%%%%%%%% Lumped Parameter 'Disc'-retization %%%%%%%


VolRad = ones([numpts_radial 1]);

InRadofslice = ones([numpts_radial numpts]);

biotslice = ones([numpts numpts_radial]);
Aradinner_array = ones(size(biotslice));
Aradouter_array = ones(size(biotslice));
darr = ones([numpts 1]);
disp('Discretizing Nozzle...')
for k = 1:length(x_metric)
    dR = (Outer_Dia(k) - Metric_Dia(k))/(2*numpts_radial);
    darr(k) = dR;
    n = 0;
    for i = 1:numpts_radial
        inner_rad_slice = Metric_Dia(k)/2 + n*dR;
        InRadofslice(i,k) = inner_rad_slice;
        Vrad = (pi)*((inner_rad_slice + dR)^2 - inner_rad_slice^2)*(Le_Metric/numpts);
        Aradinner = 2*pi*(inner_rad_slice)*(Le_Metric/numpts);
        Aradinner_array(k,i) = Aradinner;
        Aradouter = 2*pi*(inner_rad_slice + dR)*(Le_Metric/numpts);
        Aradouter_array(k,i) = Aradouter;
        VolRad(i) = Vrad;
        biotslice(k,i) = (hg_array_Nusselt(k)*Vrad)/(k_gr*Aradinner);
        n=n+1;
    end
    
end


Teatime = ones(numpts,numpts_radial,numpts_time);
disp('Iterating Heat Transfer...')
for k = 1:length(x_metric)
    hg = hg_array_Nusselt(k);
    Nmat = zeros(numpts_radial);
    Omat = zeros(numpts_radial);
    Cmat = zeros(numpts_radial);

    for i = 2:numpts_radial
        Nmat(i,i-1) = ((Aradinner_array(k,i))*k_gr)/(darr(k));
        Nmat(i,i) = -((Aradinner_array(k,i))*k_gr)/darr(k);
    end

    for i = 1:(numpts_radial-1)
        Omat(i,i+1) = (-Aradouter_array(k,i)*k_gr)/darr(k);
        Omat(i,i) = (Aradouter_array(k,i)*k_gr)/darr(k);
    end
    
    for i = 1:numpts_radial
        Cmat(i,i) = rho_gr*Cp_gr*VolRad(i);
    end
    
    Gasvec = zeros(numpts_radial,1);
    Gasvec(1) = (hg*Aradinner_array(k,1))*hotgasarray(k);

    Nmat(1,1) = -hg*Aradinner_array(k,1);
    Gmat = Nmat - Omat;
    Gmat = (inv(Cmat))*Gmat;
    Gasvec = (inv(Cmat))*Gasvec;

    for t = 1:numpts_time
        Gmatdt = (timepoints(t))*Gmat;
        Teatime(k,:,t) = (expm(Gmatdt))*(Tinit*ones([numpts_radial 1])) + (inv(Gmat))*( ...
            expm(Gmatdt) - eye(numpts_radial))*Gasvec;
    end
end


Biot_max = max(biotslice,[],"all");
disp('Maximum Biot Number:')
disp(Biot_max)

%%%%%%%%%%%%%%%%%%%%

tevalindex = numpts_time; 

Tcut = squeeze(Teatime(:,:,tevalindex))';   

% not sure how this plotting method actually works its a template but it
% looks right
Rmat = InRadofslice(:,1:numpts); % meters

[Xmat,~] = meshgrid(x_metric, 1:numpts_radial); % meters

figure(7)
surf(Xmat, Rmat, Tcut)
hold on
surf(Xmat,-Rmat,Tcut)
shading interp
view(2)
axis equal
colormap jet
xlabel('Axial position (m)')
ylabel('Radius (m)')
title(['Temperature distribution at t = ' num2str(timepoints(tevalindex)) ' s'])

%%%%%%%%%%

figure(8)

Tcut2 = squeeze(Teatime(:,1,:))';

mesh(x_metric,timepoints,Tcut2,'linewidth',1.5);
shading interp;
colormap jet;
xlabel('Axial position (m)');
ylabel('Time (s)');
zlabel('Temperature (K)');
title('Transient Nozzle Wall Temperature (K)')

%%%%%%%%%%%%%%%%%%%% Axial Position vs Temperature at end of burn %%%%%%%%%%%%%%%%%%%%

% Inner (hot-gas-side) wall temperature at the final simulated time,
% taken from the innermost radial slice (i = 1) of Tcut, which already
% corresponds to t = timepoints(tevalindex).
Twall_inner_endburn = Tcut(1,:); % K, size 1 x numpts, axial chamber -> exit

[T_peak, idx_T_peak] = max(Twall_inner_endburn);

figure(9)
plot(x_metric, Twall_inner_endburn, 'r-', 'LineWidth', 2)
hold on
plot(x_metric(idx_T_peak), T_peak, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
text(x_metric(idx_T_peak), T_peak, ...
    sprintf('  Peak = %.1f K @ x = %.4f m', T_peak, x_metric(idx_T_peak)), ...
    'FontSize',12,'VerticalAlignment','bottom')
hold off
xlabel('Axial position (m)','FontSize',14)
ylabel('Inner Wall Temperature (K)','FontSize',14)
grid on
title(['Inner Wall Temperature vs Axial Position at t = ' num2str(timepoints(tevalindex)) ' s'],'FontSize',14)

disp('Peak convective heat transfer coefficient:')
disp([num2str(hg_peak) ' W/(m^2K) at x = ' num2str(x_metric(idx_hg_peak)) ' m'])

disp('Peak inner-wall temperature at end of burn:')
disp([num2str(T_peak) ' K at x = ' num2str(x_metric(idx_T_peak)) ' m'])