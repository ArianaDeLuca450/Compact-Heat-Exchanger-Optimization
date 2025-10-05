lb1 = [0.075 0.002 1];      % lower bounds for circular (L, d_h, N)
ub1 = [0.3 0.5 6];          % upper bounds for circular
lb2 = [0.075 0.01 0.002 1]; % lower bounds for noncircular (L, a, b, N)
ub2 = [0.3 0.5 0.5 6];      % upper bounds for noncircular
intCon1 = 3;                % for circular
intCon2 = 4;                % for all others
opts = optimoptions('ga', 'Display', 'iter', 'PopulationSize', 100);
[x_opt1, fval1] = ga(@objectiveFunction1, 3, [], [], [], [], lb1, ub1, [], intCon1, opts);
[x_opt2, fval2] = ga(@objectiveFunction2, 4, [], [], [], [], lb2, ub2, [], intCon2, opts);
[x_opt3, fval3] = ga(@objectiveFunction3, 4, [], [], [], [], lb2, ub2, [], intCon2, opts);

x_opt1
x_opt2
x_opt3

function cost_circular = objectiveFunction1(x)
    % Known constants
    m = 0.023;       % mass flow rate of bleed flow (kg/s)
    p = 101325;      % pressure of the bleed flow (Pa)
    T_b = 400;       % temperature of the bleed flow (K)
    T_s = 298;       % temperture of supply air (K)
    W_target = 1465; % ideal heat transfer rate (watts)
    Pr = 0.69;       % Prandtl number
    rho = 0.8711;    % density (kg/m^3)
    cp = 1.014;      % specific heat at constant pressure (J/kgK)
    mu = 2e-5;       % dynamic viscosity (kg/ms)
    lambda = 33800;  % thermal conductivty (W/mK)
    h = 25;          % heat transfer coefficient (W/m^2K)
    delta_T = T_b - T_s; % temperature difference
    
    % Design variables
    L = x(1);        % length of duct
    d_h = x(2);      % hydraulic diameter
    N = round(x(3)); % number of ducts

    % Calculated variables for the duct
    A_c = (d_h/2)^2*pi; % flow area (m^2)
    P_w = pi*d_h;       % wetted area of pipe
    A_total = N*A_c;    % total flow area (m^2)
    A_s = P_w*L*N;      % total heat transfer area (m^2)
        

    % Calculated variables for the Flow
    v = m/(rho*A_total);      % velocity
    W = h*A_s*delta_T;        % heat transfer
    Re = (rho*v*d_h)/mu;      % Reynolds number
    St = (W*d_h)/(4*L);       % Stanton number
    Nu = 0.023*Re^0.8*Pr^0.3; % Nusselt number (Dittus-Boelter correlation)
    j = Nu/(Re*Pr^(1/3));     % Colburn factor
    f = 16/Re;                % fanning friction factor
    delta_p = f*(L/d_h)*((rho*v^2)/2);    % pressure drop
    G = sqrt(((St/f)/W)*(2*rho*delta_p)); % core mass velocity (kg/m^2s)
    cost_circular = (abs(W_target - W))/W_target + d_h*0.5 + delta_p*2e-3;
end

function cost_sine = objectiveFunction2(x)
    % Known constants
    m = 0.023;       % mass flow rate of bleed flow (kg/s)
    p = 101325;      % pressure of the bleed flow (Pa)
    T_b = 400;       % temperature of the bleed flow (K)
    T_s = 298;       % temperture of supply air (K)
    W_target = 1465; % ideal heat transfer rate (watts)
    Pr = 0.69;       % Prandtl number
    rho = 0.8711;    % density (kg/m^3)
    cp = 1.014;      % specific heat at constant pressure (J/kgK)
    mu = 2e-5;       % dynamic viscosity (kg/ms)
    lambda = 33800;  % thermal conductivty (W/mK)
    h = 25;          % heat transfer coefficient (W/m^2K)
    delta_T = T_b - T_s; % temperature difference
    
    % Design variables
    L = x(1);        % length of duct
    a = x(2);        % spread
    b = x(3);        % mean
    N = round(x(4)); % number of ducts

    % Calculated variables for the duct
    % flow area (m^2)
    A_c = integral(@(x) b*(1 + cos(pi*x/a)), -a, a);
    % wetted area of pipe
    P_w = integral(@(x) sqrt(1 + (b*pi/a)^2 * sin(pi*x/a).^2), -a, a);
    A_total = N*A_c;    % total flow area (m^2)
    A_s = P_w*L*N;      % total heat transfer area (m^2)
    d_h = (4*A_c)/P_w;  % hydraulic diameter
        

    % Calculated variables for the Flow
    v = m/(rho*A_total);      % velocity
    W = h*A_s*delta_T;        % heat transfer
    Re = (rho*v*d_h)/mu;      % Reynolds number
    St = (W*d_h)/(4*L);       % Stanton number
    Nu = 0.023*Re^0.8*Pr^0.3; % Nusselt number (Dittus-Boelter correlation)
    j = Nu/(Re*Pr^(1/3));     % Colburn factor
    f = 16/Re;                % fanning friction factor
    delta_p = f*(L/d_h)*((rho*v^2)/2);    % pressure drop
    G = sqrt(((St/f)/W)*(2*rho*delta_p)); % core mass velocity (kg/m^2s)
    cost_sine = (abs(W_target - W))/W_target + d_h*0.5 + delta_p*2e-3;
end

function cost_rect = objectiveFunction3(x)
    % Known constants
    m = 0.023;       % mass flow rate of bleed flow (kg/s)
    p = 101325;      % pressure of the bleed flow (Pa)
    T_b = 400;       % temperature of the bleed flow (K)
    T_s = 298;       % temperture of supply air (K)
    W_target = 1465; % ideal heat transfer rate (watts)
    Pr = 0.69;       % Prandtl number
    rho = 0.8711;    % density (kg/m^3)
    cp = 1.014;      % specific heat at constant pressure (J/kgK)
    mu = 2e-5;       % dynamic viscosity (kg/ms)
    lambda = 33800;  % thermal conductivty (W/mK)
    h = 25;          % heat transfer coefficient (W/m^2K)
    delta_T = T_b - T_s; % temperature difference
    
    % Design variables
    L = x(1);        % length of duct
    a = x(2);        % height
    b = x(3);        % width
    N = round(x(4)); % number of ducts

    % Calculated variables for the duct 
    A_c = a*b; % flow area (m^2)
    d_h = (4*a*b)/(a + b);  % hydraulic diameter
    P_w = pi*d_h;           % wetted area of pipe
    A_total = N*A_c;        % total flow area (m^2)
    A_s = P_w*L*N;          % total heat transfer area (m^2)
        

    % Calculated variables for the Flow
    v = m/(rho*A_total);      % velocity
    W = h*A_s*delta_T;        % heat transfer
    Re = (rho*v*d_h)/mu;      % Reynolds number
    St = (W*d_h)/(4*L);       % Stanton number
    Nu = 0.023*Re^0.8*Pr^0.3; % Nusselt number (Dittus-Boelter correlation)
    j = Nu/(Re*Pr^(1/3));     % Colburn factor
    f = 16/Re;                % fanning friction factor
    delta_p = f*(L/d_h)*((rho*v^2)/2);    % pressure drop
    G = sqrt(((St/f)/W)*(2*rho*delta_p)); % core mass velocity (kg/m^2s)
    cost_rect = (abs(W_target - W))/W_target + d_h*0.5 + delta_p*2e-3;
end