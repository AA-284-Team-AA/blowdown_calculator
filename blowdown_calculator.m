%% Tank Blowdown Analysis

% Thomas White 2/23/2020

% Outputs:
% An idea of the flow rate over time. 
% Initial and final pressure in the tank. 


%% Unit Conversion

% Unit Conversion
psi_to_Pa = 6894.75729; % 1 psi in Pa
in_to_m = 0.0254; % 1 in in m
mm_to_m = 1e-3; % 1 mm in m
lbf_to_N = 4.44822162; % 1 lbf in N
lbm_to_kg = 0.453592; % 1 lbm in kg
atm_to_Pa = 101325; % 1 atm in Pa
L_to_m3 = 1e-3; % 1 L in m^3


%% Inputs

% Assumptions - flow is choked at the regulator. 

gamma = 1.4;
R = 8.31; % Universal gas constant
Mw = 0.0289; % kg/mol Molar weight of gas
T0 = 298; %K initial 
P0 = 1000*psi_to_Pa; % Pa
V_tank = 0.026; % m^3 
% Cd = 1;
% A = 2e-6; % m^2
Cv = 0.06

t = linspace(0, 20);


%% Calc Flow Velocity Out of Tank

% Choked in regulator at M = 1

c = sqrt(((2*gamma)/(gamma + 1))*(R*T0)/Mw)


% Find flow coefficient Cd*A

CdA = Cv * 16.2*mm_to_m^2;


%% Create Dischange Coefficient 

c0 = sqrt(gamma*R*T0/Mw)

tau = (V_tank/(CdA*c0))*((gamma + 1)/2)^((gamma + 1)/(2*(gamma-1)))


%% Find results

% Isothermal Tank

P_tank_iso = P0*exp(-t/tau);


figure 
hold on
plot(t, P_tank_iso/psi_to_Pa)

% Adiabatic Tank

P_tank_adi = P0*(1 + ((gamma - 1)/2)*(-t/tau)).^((2*gamma)/(gamma - 1));

plot(t, P_tank_adi/psi_to_Pa)

title("Run Tank Blowdown")
legend("Isothermal Blowdown", "Adiabatic Blowdown")
xlabel("Time (s)")
ylabel("Tank Pressure (psi)")
