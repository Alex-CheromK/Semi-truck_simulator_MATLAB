% This MATLAB script uses symbolic processing to analytically find the
% equations of motion describing the dynamics of a truck-trailer system
% with non-holonomically constrained motion.

%% Setup workspace
clear;
clc;

%% Define symbolic variables
% Truck local velocities and angular velocity
syms vAx vAy omegaA real

% Trailer angle relative to truck
syms omegaB gamma real

% Truck local forward acceleration
syms vAx_dot real

% Truck front wheel angle and angular velocity
syms phi phi_dot real

% Truck thrust force
syms T real

% System lengths...
% L1 - From truck CG to front wheels
% L2 - From truck CG to rear wheels/pivot point
% L3 - From pivot point to trailer CG
% L4 - From pivot point to trailer rear wheels
syms L1 L2 L3 L4 real

% Lateral wheel forces...
% lambda1 - Truck front wheels
% lambda2 - Truck rear wheels
% lambda3 - Trailer rear wheels
syms lambda1 lambda2 lambda3 real

% Truck and trailer masses and moments of inertia
syms mA mB IxxA IyyA IzzA IxxB IyyB IzzB real

% Global properties...
% g - Gravitational acceleration constant
% mu_k - Coefficient of kinetic friction between tires and road
% rho - Air density
syms g mu_k rho real

% Vehicle drag coefficient and drag refernce area
syms C_D A real

%% Kinematics
% System position vectors
r_A_1 = [L1 0 0]';
r_A_2 = [-L2 0 0]';
r_2_3 = [-L4*cos(gamma) -L4*sin(gamma) 0]';
r_2_B = [-L3*cos(gamma) -L3*sin(gamma) 0]';
r_2_1 = [(L1 + L2) 0 0]';
r_2_A = -r_A_2;

% Moment of inertia matrices
I_A = [
    IxxA 0 0;
    0 IyyA 0;
    0 0 IzzA];
I_B = [
    IxxB 0 0;
    0 IyyB 0;
    0 0 IzzB];

% Angular velocity and velocity vectors of truch and trailer CGs
omegaA_vec = [0 0 omegaA]';
omegaB_vec = [0 0 omegaB]';
vA_vec = [vAx vAy 0]';
vB_vec =...
    vA_vec +...
    cross(omegaB_vec,r_2_B) +...
    cross(omegaA_vec,r_A_2 + r_2_B);

% Velocity vectors of other points
v1_vec = vA_vec + cross(omegaA_vec,r_A_1);
v2_vec = vA_vec + cross(omegaA_vec,r_A_2);
v3_vec =...
    vA_vec +...
    cross(omegaB_vec,r_2_3) +...
    cross(omegaA_vec,r_A_2 + r_2_3);

% Unit vector pointing in the directions of lateral wheel forces
u1_vec = [-sin(phi) cos(phi) 0]';
u2_vec = [0 1 0]';
u3_vec = [-sin(gamma) cos(gamma) 0]';

% Define a system of equations that represent the non-holonomic constraints
Eqs = [
    dot(v1_vec,u1_vec) == 0;
    dot(v2_vec,u2_vec) == 0;
    dot(v3_vec,u3_vec) == 0];
% Solve the system of equations for...
Sol = solve(Eqs,[vAy,omegaA,omegaB]);
% ...truck lateral velocity...
vAy = simplify(Sol.vAy);
% ...truck angular velocity...
omegaA = simplify(Sol.omegaA);
% ...and trailer angular velocity...
omegaB = simplify(Sol.omegaB);
% ... all as functions of vAx (truck forward speed)

% Differentiate above expressions
vAy_dot = (L2/(L1 + L2))*(vAx_dot*tan(phi) + vAx*phi_dot*sec(phi)^2);
omegaA_dot = (1/(L1 + L2))*(vAx_dot*tan(phi) + vAx*phi_dot*sec(phi)^2);
omegaB_dot = (-1/(L4*(L1 + L2)))*...
    ((L4*phi_dot*sec(phi)^2 + (L1 + L2)*omegaB*cos(gamma))*vAx +...
    (L4*tan(phi) + (L1 + L2)*sin(gamma))*vAx_dot);

% Compile acceleration vectors for EOMs
vA_dot_vec = [vAx_dot vAy_dot 0]';
alphaA_vec = [0 0 omegaA_dot]';
alphaB_vec = [0 0 omegaB_dot]';
aA_vec = vA_dot_vec + cross(omegaA_vec,vA_vec);
aB_vec = simplify(...
    aA_vec +...
    cross(alphaB_vec,r_2_B) +...
    cross(omegaB_vec,cross(omegaB_vec,r_2_B) + cross(omegaA_vec,r_2_B)) +...
    cross(alphaA_vec,r_A_2 + r_2_B) +...
    cross(omegaA_vec,cross(omegaB_vec,r_2_B) + cross(omegaA_vec,r_A_2 + r_2_B)));

% Define normal forces at each wheel set for friction calculations
N1 = L2/(L1 + L2)*mA*g;
N2 = L1/(L1 + L2)*mA*g + (L4 - L3)/L4*mB*g;
N3 = L3/L4*mB*g;

% Define external forces
syms F_Dx F_Dy Ff1x Ff1y Ff2x Ff2y Ff3x Ff3y real
F_D = [F_Dx F_Dy 0]';
Ff1 = [Ff1x Ff1y 0]';
Ff2 = [Ff2x Ff2y 0]';
Ff3 = [Ff3x Ff3y 0]';
% Note: Above forces are simplified to components for now since these
% components do not contain any dependent variables and can easily be
% replaced with numerical values later on. There is no need to
% unnecessarily complicate the symbolic processing.
F_T = T*[cos(phi) sin(phi) 0]';
R1 = [-lambda1*sin(phi) lambda1*cos(phi) 0]';
R2 = [0 lambda2 0]';
R3 = [-lambda3*sin(gamma) lambda3*cos(gamma) 0]';

% Define rates of change of linear...
dL_A_dt = mA*aA_vec;
dL_B_dt = mB*aB_vec;
% ...and angular momenta
dH_A_dt = I_A*alphaA_vec + cross(omegaA_vec,I_A*omegaA_vec);
dH_B_dt = I_B*alphaB_vec + cross(omegaB_vec,I_B*omegaB_vec);

% Define the EOMs as a system of equations
Eqs = [
    F_T + F_D + Ff1 + Ff2 + Ff3 + R1 + R2 + R3 == dL_A_dt + dL_B_dt;
    cross(r_2_1,F_T + Ff1 + R1) - cross(r_2_A,dL_A_dt) == dH_A_dt;
    cross(r_2_3,R3 + Ff3) - cross(r_2_B,dL_B_dt) == dH_B_dt];
Eqs = Eqs([1 2 6 9]);
% Solve above system for...
Sol = solve(Eqs,[vAx_dot,lambda1,lambda2,lambda3]);
% ...truch forward acceleration...
vAx_dot = simplify(Sol.vAx_dot);
% ...and wheel lateral forces
lambda1 = simplify(Sol.lambda1);
lambda2 = simplify(Sol.lambda2);
lambda3 = simplify(Sol.lambda3);