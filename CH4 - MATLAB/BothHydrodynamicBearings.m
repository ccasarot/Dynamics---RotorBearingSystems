%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project 2 - Chapter 4                                 %
%                                                       %
%               Copenhagen, Spring semester 2023        %
%                                                       %
%                     Christian Casarotto - s223302     %
%                                                       %
% Using now both hydrodynamic bearings                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONTENT
% In this file the system's Campbell's diagram and stability map is studied
% again for the configuration with both hydrodynamic bearings.

close all
clear all

% Bearing Properties 
% Table 1a : Two-axial-groove bearing, L/D = 0.5
%
%        S     E     Phi   Q     P     T   Kxx   Kxy   Kyx   Kyy  Bxx   Bxy  Byx  Byy
%
Table=[6.430 0.071 81.89 0.121 0.860  5.7  1.55 14.41 -6.60 1.88 28.75 1.89 1.89 13.31
       3.937 0.114 77.32 0.192 0.846  5.9  1.57  9.27 -4.20 1.89 18.44 1.93 1.93  8.58
       2.634 0.165 72.36 0.271 0.833  6.2  1.61  6.74 -3.01 1.91 13.36 2.00 2.00  6.28
       2.030 0.207 68.75 0.332 0.835  6.6  1.65  5.67 -2.50 1.93 11.18 2.07 2.07  5.33
       1.656 0.244 65.85 0.383 0.835  7.0  1.69  5.06 -2.20 1.95  9.93 2.15 2.15  4.80
       0.917 0.372 57.45 0.540 0.850  8.5  2.12  4.01 -1.30 1.85  7.70 2.06 2.06  3.23     
       0.580 0.477 51.01 0.651 0.900 10.5  2.67  3.70 -0.78 1.75  6.96 1.94 1.94  2.40 
       0.378 0.570 45.43 0.737 0.977 13.4  3.33  3.64 -0.43 1.68  6.76 1.87 1.87  1.89
       0.244 0.655 40.25 0.804 1.096 17.9  4.21  3.74 -0.13 1.64  6.87 1.82 1.82  1.54
       0.194 0.695 37.72 0.833 1.156 21.3  4.78  3.84  0.01 1.62  7.03 1.80 1.80  1.40
       0.151 0.734 35.20 0.858 1.240 25.8  5.48  3.98  0.15 1.61  7.26 1.79 1.79  1.27
       0.133 0.753 33.93 0.870 1.289 28.7  5.89  4.07  0.22 1.60  7.41 1.79 1.79  1.20
       0.126 0.761 33.42 0.875 1.310 30.0  6.07  4.11  0.25 1.60  7.48 1.79 1.79  1.18
       0.116 0.772 32.65 0.881 1.343 32.2  6.36  4.17  0.30 1.60  7.59 1.79 1.79  1.15
       0.086 0.809 30.04 0.902 1.473 41.4  7.51  4.42  0.47 1.59  8.03 1.79 1.79  1.03
       0.042 0.879 24.41 0.936 1.881 80.9 11.45  5.23  0.92 1.60  9.48 1.80 1.80  0.82 ];
 
S_values = Table(:, 1);
E_values = Table(:, 2);
Phi_values = Table(:, 3);
Q_values = Table(:, 4);
P_values = Table(:, 5);
T_values = Table(:, 6);
Kxx_values = Table(:, 7);     Bxx_values = Table(:, 11);
Kxy_values = Table(:, 8);     Bxy_values = Table(:, 12);
Kyx_values = Table(:, 9);     Byx_values = Table(:, 13);
Kyy_values = Table(:, 10);    Byy_values = Table(:, 14);

N_campbell=500;

for iii=1:N_campbell
    
  %Omega= (150/60/300*iii*2*pi)+((5100*2*pi)/60); % angular velocity [rad/s]
  Omega = pi+(iii)*2*pi;                    % angular velocity [rad/s]
  Omegarpm(iii) = (Omega*60)/(2*pi);   % angular velocity [rpm]
  N=Omega/(2*pi);                      % Ang vel in Hz
  N_vector(iii) = N;  % collect for plotting
       


%% FLUID BEARING 1 (sx) % % % % % % % % % % % % % % % % % % % % % % % % % % 



% Definition of Nondimensional Bearing Parameters
W = 892.7;                        % W - external load in N; 
eta = 0.0277 * exp(1)^(0.034*(40-55)); % viscosity % η - oil viscosity

% Cross section

C = 1e-4;                              % [m] clearance % C is the bearing clearance in m; 
d = 99.6/1000;                         % [m] diameter of the shaft % D is the bearing inner diameter in m;
r = d/2; 
R = C + r; % form C = R - r
D = R*2;   % from R = D/2;
syms L
L = solve(L/D==0.5,L);                 % as L/d is 0.5 % L is the bearing width in m; 

S = double(eta*N*L*D/W*(R/C)^2);               % S - Sommerfeld Number

% Perform spline interpolation
E_interpolated = interp1(S_values, E_values, S, 'spline');
Phi_interpolated = interp1(S_values, Phi_values, S, 'spline');
Q_interpolated = interp1(S_values, Q_values, S, 'spline');
P_interpolated = interp1(S_values, P_values, S, 'spline');
T_interpolated = interp1(S_values, T_values, S, 'spline');
Kxx_interpolated = interp1(S_values, Kxx_values, S, 'spline');
Kxy_interpolated = interp1(S_values, Kxy_values, S, 'spline');
Kyx_interpolated = interp1(S_values, Kyx_values, S, 'spline');
Kyy_interpolated = interp1(S_values, Kyy_values, S, 'spline');
Bxx_interpolated = interp1(S_values, Bxx_values, S, 'spline');
Bxy_interpolated = interp1(S_values, Bxy_values, S, 'spline');
Byx_interpolated = interp1(S_values, Byx_values, S, 'spline');
Byy_interpolated = interp1(S_values, Byy_values, S, 'spline');

% Coefficients of the Stiffness Matrix  
kxx_1 =  Kxx_interpolated*W/C  ; % [N/m]  
kxy_1 =  Kxy_interpolated*W/C  ; % [N/m]
kyx_1 =  Kyx_interpolated*W/C  ; % [N/m]
kyy_1 =  Kyy_interpolated*W/C ; % [N/m]

% Coefficients of the Damping Matrix 
dxx_1 =  Bxx_interpolated*W/(Omega*C)  ; % [N/(m/s)]
dxy_1 =  Bxy_interpolated*W/(Omega*C)  ; % [N/(m/s)]
dyx_1 =  Byx_interpolated*W/(Omega*C)  ; % [N/(m/s)]
dyy_1 =  Byy_interpolated*W/(Omega*C)  ; % [N/(m/s)]



%% FLUID BEARING 2 (sx) % % % % % % % % % % % % % % % % % % % % % % % % % % 



% Definition of Nondimensional Bearing Parameters
W = 45;                        % W - external load in N; 
eta = 0.0277 * exp(1)^(0.034*(40-55)); % viscosity % η - oil viscosity

% Cross section
C = 0.85e-4;                           % [m] clearance % C is the bearing clearance in m; 
d = 50/1000;                         % [m] diameter of the shaft % D is the bearing inner diameter in m;
r = d/2; 
R = C + r; % form C = R - r
D = R*2;   % from R = D/2;
syms L
L = solve(L/D==0.5,L);                 % as L/d is 0.5 % L is the bearing width in m; 

S = double(eta*N*L*D/W*(R/C)^2);               % S - Sommerfeld Number

% Perform spline interpolation
E_interpolated = interp1(S_values, E_values, S, 'spline');
Phi_interpolated = interp1(S_values, Phi_values, S, 'spline');
Q_interpolated = interp1(S_values, Q_values, S, 'spline');
P_interpolated = interp1(S_values, P_values, S, 'spline');
T_interpolated = interp1(S_values, T_values, S, 'spline');
Kxx_interpolated = interp1(S_values, Kxx_values, S, 'spline');
Kxy_interpolated = interp1(S_values, Kxy_values, S, 'spline');
Kyx_interpolated = interp1(S_values, Kyx_values, S, 'spline');
Kyy_interpolated = interp1(S_values, Kyy_values, S, 'spline');
Bxx_interpolated = interp1(S_values, Bxx_values, S, 'spline');
Bxy_interpolated = interp1(S_values, Bxy_values, S, 'spline');
Byx_interpolated = interp1(S_values, Byx_values, S, 'spline');
Byy_interpolated = interp1(S_values, Byy_values, S, 'spline');

% Coefficients of the Stiffness Matrix  
kxx_2 =  Kxx_interpolated*W/C  ; % [N/m]  
kxy_2 =  Kxy_interpolated*W/C  ; % [N/m]
kyx_2 =  Kyx_interpolated*W/C  ; % [N/m]
kyy_2 =  Kyy_interpolated*W/C ; % [N/m]

% Coefficients of the Damping Matrix 
dxx_2 =  Bxx_interpolated*W/(Omega*C)  ; % [N/(m/s)]
dxy_2 =  Bxy_interpolated*W/(Omega*C)  ; % [N/(m/s)]
dyx_2 =  Byx_interpolated*W/(Omega*C)  ; % [N/(m/s)]
dyy_2 =  Byy_interpolated*W/(Omega*C)  ; % [N/(m/s)]
    


%% MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



 RedFactorId = 0.9;
 E    = 2.0E11;      % {elasticity modulus [N/m^2}
 RAco = 7800;        % {steel density [kg/m^3]}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DEFINITION OF THE STRUCTURE OF THE MODEL   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 NE=23;         % number of shaft elements
 GL = (NE+1)*4; % number of degree of freedom
 ND=2;          % number of discs
 NM=2;          % number of bearings
 CD1=4;         % node - disc 1
 CD2=6;         % node - disc 2
 CMM1=9;        % node - bearing 1
 CMM2=NE-1;     % node - bearing 2
 
% (A) DISK 2
 Rd_2   = 295/2000;          % external radius of the disc [m]
 Ri_2   = 70/2000;        % internal radius of the disc [m]
 espD_2 = 100/1000 ;         % disc thickness  [m]
 % Removed ring
 Rr_2 = 271/2000;             % Internal ratius ring
 Lr_2 = (35/3)/1000;          % Lenght of the ring
 Rdist_2 = 35/1000;           % Ring distance
 % Other rings
 R_external_SmallRing_2 = 90/2000;             % Internal ratius ring
 L_SmallRing_2 = 20/1000;          % Lenght of the ring
 Rdist_SmallRing_2 = 15/1000;           % Ring distance
 R_external_BigRing_2 = 110/2000;             % Internal ratius ring
 L_BigRing_2 = 25/1000;          % Lenght of the ring
 Rdist_BigRing_2 = 37.5/1000;           % Ring distance

 % disc mass [kg]  
     Mde_2 = pi*Rd_2^2*espD_2*RAco;  % Mass disc external
     Mdi_2 = pi*Ri_2^2*espD_2*RAco;  % Mass disc internal
     Mre_2 = pi*Rd_2^2*Lr_2*RAco;    % Mass ring external
     Mri_2 = pi*Rr_2^2*Lr_2*RAco;    % Mass ring internal
 M_external_SmallRing_2 = pi*R_external_SmallRing_2^2*L_SmallRing_2*RAco;    
 M_internal_SmallRing_2 = pi*Ri_2^2*L_SmallRing_2*RAco;    
 M_external_BigRing_2 = pi*R_external_BigRing_2^2*L_BigRing_2*RAco;    
 M_internal_BigRing_2 = pi*Ri_2^2*L_BigRing_2*RAco;     
     MasD_2 = Mde_2 - Mdi_2 - 3*(Mre_2 - Mri_2) - ...
              (M_external_SmallRing_2 - M_internal_SmallRing_2) - ...
              (M_external_BigRing_2 - M_internal_BigRing_2);

 % Inertia [Kgm^2]
 Id_2 = ((1/4)*Rd_2^2 + (1/12)*espD_2^2) * Mde_2 ... % Overall ext
       -((1/4)*Ri_2^2 + (1/12)*espD_2^2) * Mdi_2 ... % Overall int
           -((1/4)*Rd_2^2 + (1/12)*Lr_2^2) * Mre_2 ...   % Central ring ext
           +((1/4)*Rr_2^2 + (1/12)*Lr_2^2) * Mri_2 ...   % Central ring int
       - 2 * ((1/4)*Rd_2^2 + (1/12)*Lr_2^2) * Mre_2 ...   % Side rings ext
       + 2 * ((1/4)*Rr_2^2 + (1/12)*Lr_2^2) * Mri_2 ...   % Side rings int
       - 2 * (Mre_2-Mri_2) * Rdist_2^2 ...                   % Transport 
           - ((1/4)*R_external_SmallRing_2^2 + (1/12)*L_SmallRing_2^2) * M_external_SmallRing_2 ...   
           + ((1/4)*Ri_2^2 + (1/12)*L_SmallRing_2^2) * M_internal_SmallRing_2 ...   
           - (M_external_SmallRing_2-M_internal_SmallRing_2) * Rdist_SmallRing_2^2 ...
           - ((1/4)*R_external_BigRing_2^2 + (1/12)*L_BigRing_2^2) * M_external_BigRing_2 ...   
           + ((1/4)*Ri_2^2 + (1/12)*L_BigRing_2^2) * M_internal_BigRing_2 ...   
           - (M_external_BigRing_2-M_internal_BigRing_2) * Rdist_BigRing_2^2;

 Ip_2 =   1/2*(Rd_2^2)*Mde_2 ... % Overall disk
        - 1/2*(Ri_2^2)*Mdi_2 ...
            - 1/2*(Rd_2^2)*Mre_2 ... % External ring (scanalature)
            + 1/2*(Rr_2^2)*Mri_2 ...
        - 1/2*(R_external_SmallRing_2^2)*M_external_SmallRing_2 ... % Internal ring piccolo (sedi)
        + 1/2*(Ri_2^2)*M_internal_SmallRing_2 ...
            - 1/2*(R_external_BigRing_2^2)*M_external_BigRing_2 ... % Internal ring grande (sedi)
            + 1/2*(Ri_2^2)*M_internal_BigRing_2;

 % (C) SHAFT  
 Rint = 0;             % shaft internal radius [m]                      

 % shaft elements [m]
 l(1)  = 25/1000;      rx(1)  = 30/2000;
 l(2)  = 40/1000;      rx(2)  = 65/2000;
 l(3)  = 35/1000;      rx(3)  = 70/2000;
 l(4)  = 40/1000;      rx(4)  = 70/2000;
 l(5)  = 50/1000;      rx(5)  = 70/2000;
 l(6)  = 50/1000;      rx(6)  = 70/2000;
 l(7)  = 55/1000;      rx(7)  = 90/2000;
 l(8)  = 60/1000;      rx(8)  = 99.6/2000;
 l(9)  = 60/1000;      rx(9)  = 99.6/2000;
 l(10) = 62.5/1000;    rx(10) = 90/2000;
 l(11) = 62.5/1000;    rx(11) = 90/2000;
 l(12) = 62.5/1000;    rx(12) = 90/2000;
 l(13) = 62.5/1000;    rx(13) = 90/2000;
 l(14) = 62.5/1000;    rx(14) = 90/2000;
 l(15) = 62.5/1000;    rx(15) = 90/2000;
 l(16) = 62.5/1000;    rx(16) = 90/2000;
 l(17) = 62.5/1000;    rx(17) = 90/2000;
 l(18) = 62.5/1000;    rx(18) = 90/2000;
 l(19) = 62.5/1000;    rx(19) = 90/2000;
 l(20) = 35/1000;      rx(20) = 70/2000; 
 l(21) = 11.5/1000;    rx(21) = 50/2000; 
 l(22) = 11.5/1000;    rx(22) = 50/2000; 
 l(23) = 52/1000;      rx(23) = 40/2000; 

% internal radius of shaft elements [m]
    for i=1:NE                 ri(i)=Rint; end
% density of shaft elements [kg/m]
    for i=1:NE                 ro(i) = RAco; end
% transversal areal of the shaft elements [m^2]}
    for i=1:NE St(i) = pi*(rx(i)^2-ri(i)^2); end
% area moment of inertia of the shaft elements [m^4]}
    for i=1:NE II(i)=pi*(rx(i)^4-ri(i)^4)/4; end
    
% MOUNTING THE GLOBAL MATRICES          
M=zeros(GL); G=zeros(GL); K=zeros(GL); Damp=zeros(GL);   
   
% Mass matrices of shaft elements due to linear and angular movements
a=1; b=8;
for n=1:NE 

  MteAux= [156       0         0          22*l(n)    54        0         0          -13*l(n)
           0         156       -22*l(n)   0          0         54        13*l(n)    0
           0         -22*l(n)  4*l(n)^2   0          0         -13*l(n)  -3*l(n)^2  0
           22*l(n)   0         0          4*l(n)^2   13*l(n)   0         0          -3*l(n)^2
           54        0         0          13*l(n)    156       0         0          -22*l(n)
           0         54        -13*l(n)   0          0         156       22*l(n)    0
           0         13*l(n)   -3*l(n)^2  0          0         22*l(n)   4*l(n)^2   0
           -13*l(n)  0         0          -3*l(n)^2  -22*l(n)  0         0          4*l(n)^2];

  Mte = ((ro(n)*St(n)*l(n))/420)*MteAux;
  
  MreAux= [36      0        0         3*l(n)    -36      0       0         3*l(n)
           0       36       -3*l(n)   0         0        -36     -3*l(n)   0
           0       -3*l(n)  4*l(n)^2  0         0        3*l(n)  -l(n)^2   0
           3*l(n)  0        0         4*l(n)^2  -3*l(n)  0       0         -l(n)^2
           -36     0        0         -3*l(n)   36       0       0         -3*l(n)
           0       -36      3*l(n)    0         0        36      3*l(n)    0
           0       -3*l(n)  -l(n)^2   0         0        3*l(n)  4*l(n)^2  0
           3*l(n)  0        0         -l(n)^2   -3*l(n)  0       0         4*l(n)^2];

  Mre = ((ro(n)*St(n)*(rx(n)^2-ri(n)^2))/(120*l(n)))*MreAux;
 
  MauxT=Mte+Mre;

   for f=a:b
    for g=a:b
     M(f,g)=M(f,g)+MauxT(f-(n-1)*4,g-(n-1)*4);
    end
   end
a=a+4; b=b+4;
end

% Adding the mass matrices of DISC 2 - Only disc 2
   
   M((CD2-1)*4+1,(CD2-1)*4+1)=M((CD2-1)*4+1,(CD2-1)*4+1)+MasD_2;
   M((CD2-1)*4+2,(CD2-1)*4+2)=M((CD2-1)*4+2,(CD2-1)*4+2)+MasD_2;
   M((CD2-1)*4+3,(CD2-1)*4+3)=M((CD2-1)*4+3,(CD2-1)*4+3)+Id_2*RedFactorId;
   M((CD2-1)*4+4,(CD2-1)*4+4)=M((CD2-1)*4+4,(CD2-1)*4+4)+Id_2*RedFactorId;
               
% Gyroscopic matrix of shaft elements
a=1; b=8;
for n=1:NE

   GeAux=[0        -36       3*l(n)   0          0        36       3*l(n)    0
          36       0        0         3*l(n)     -36      0        0         3*l(n)
          -3*l(n)  0        0         -4*l(n)^2  3*l(n)   0        0         l(n)^2
          0        -3*l(n)  4*l(n)^2  0          0        3*l(n)   -l(n)^2   0
          0        36       -3*l(n)   0          0        -36      -3*l(n)   0
          -36      0        0         -3*l(n)    36       0        0         -3*l(n)
          -3*l(n)  0        0         l(n)^2     3*l(n)   0        0         -4*l(n)^2
          0        -3*l(n)  -l(n)^2   0          0        3*l(n)   4*l(n)^2  0        ];

   Ge = 2*((ro(n)*St(n)*(rx(n)^2+ri(n)^2))/(120*l(n)))*GeAux;
   
   for f=a:b
    for g=a:b
     G(f,g)=G(f,g)+Ge(f-(n-1)*4,g-(n-1)*4);
    end
   end
a=a+4; b=b+4;
end

% Adding the gyroscopic matrices of DISC 2 - Only disc 2
   
   G((CD2-1)*4+3,(CD2-1)*4+4)=G((CD2-1)*4+3,(CD2-1)*4+4)-Ip_2;
   G((CD2-1)*4+4,(CD2-1)*4+3)=G((CD2-1)*4+4,(CD2-1)*4+3)+Ip_2;
           
% Stiffness matrix of shaft elements due to bending
a=1; b=8;
for n=1:NE

  KbeAux= [12      0        0         6*l(n)    -12      0       0         6*l(n)
           0       12       -6*l(n)   0         0        -12     -6*l(n)   0
           0       -6*l(n)  4*l(n)^2  0         0        6*l(n)  2*l(n)^2  0
           6*l(n)  0        0         4*l(n)^2  -6*l(n)  0       0         2*l(n)^2
           -12     0        0         -6*l(n)   12       0       0         -6*l(n)
           0       -12      6*l(n)    0         0        12      6*l(n)    0
           0       -6*l(n)  2*l(n)^2  0         0        6*l(n)  4*l(n)^2  0
           6*l(n)  0        0         2*l(n)^2  -6*l(n)  0       0         4*l(n)^2];

  Kbe = ((E*II(n))/(l(n)^3))*KbeAux;

  for f=a:b
   for g=a:b
    K(f,g)=K(f,g)+Kbe(f-(n-1)*4,g-(n-1)*4);
   end
  end
a=a+4; b=b+4;
end

% Adding the stiffness on the bearing position

   % Hydro bearing 1
   K((CMM1-1)*4+1,(CMM1-1)*4+1)=K((CMM1-1)*4+1,(CMM1-1)*4+1) + kxx_1;
   K((CMM1-1)*4+1,(CMM1-1)*4+2)=K((CMM1-1)*4+1,(CMM1-1)*4+2) + kxy_1;
   K((CMM1-1)*4+2,(CMM1-1)*4+1)=K((CMM1-1)*4+2,(CMM1-1)*4+1) + kyx_1;
   K((CMM1-1)*4+2,(CMM1-1)*4+2)=K((CMM1-1)*4+2,(CMM1-1)*4+2) + kyy_1;
   
   % Hydro bearing 2
   K((CMM2-1)*4+1,(CMM2-1)*4+1)=K((CMM2-1)*4+1,(CMM2-1)*4+1) + kxx_2;
   K((CMM2-1)*4+1,(CMM2-1)*4+2)=K((CMM2-1)*4+1,(CMM2-1)*4+2) + kxy_2;
   K((CMM2-1)*4+2,(CMM2-1)*4+1)=K((CMM2-1)*4+2,(CMM2-1)*4+1) + kyx_2;
   K((CMM2-1)*4+2,(CMM2-1)*4+2)=K((CMM2-1)*4+2,(CMM2-1)*4+2) + kyy_2;

% Adding the damping on the bearing position

   Damp((CMM1-1)*4+1,(CMM1-1)*4+1)=Damp((CMM1-1)*4+1,(CMM1-1)*4+1) + dxx_1;
   Damp((CMM1-1)*4+1,(CMM1-1)*4+2)=Damp((CMM1-1)*4+1,(CMM1-1)*4+2) + dxy_1;
   Damp((CMM1-1)*4+2,(CMM1-1)*4+1)=Damp((CMM1-1)*4+2,(CMM1-1)*4+1) + dyx_1;
   Damp((CMM1-1)*4+2,(CMM1-1)*4+2)=Damp((CMM1-1)*4+2,(CMM1-1)*4+2) + dyy_1;

   Damp((CMM2-1)*4+1,(CMM2-1)*4+1)=Damp((CMM2-1)*4+1,(CMM2-1)*4+1) + dxx_2;
   Damp((CMM2-1)*4+1,(CMM2-1)*4+2)=Damp((CMM2-1)*4+1,(CMM2-1)*4+2) + dxy_2;
   Damp((CMM2-1)*4+2,(CMM2-1)*4+1)=Damp((CMM2-1)*4+2,(CMM2-1)*4+1) + dyx_2;
   Damp((CMM2-1)*4+2,(CMM2-1)*4+2)=Damp((CMM2-1)*4+2,(CMM2-1)*4+2) + dyy_2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    GLOBAL MATHEMATICAL MODEL                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mglob = [        M        zeros(size(M))   ; % form the other file of ch4
            zeros(size(M))            M  ] ;
Kglob = [   -G*Omega+Damp         K         ; 
           -M                zeros(size(M))];

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %              MODAL ANALYSIS                  %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Calculating Eigenvectors and Eigenvalues

 [U,lambda]=eig(-Kglob,Mglob);
 lambda_d=diag(lambda);
 [lam,p]=sort(diag(abs(imag(lambda))));
 %U=U(:,p); % nothing changes?
 
 nnn=192;
 
     lambda_d(p(1:nnn))/(2*pi);         
     lambda_campbell(iii,:)=imag(diag(lambda(p(1:nnn),p(1:nnn))))'/2/pi;
     lambda_stability(iii,:)=real(diag(lambda(p(1:nnn),p(1:nnn))))';  
    
 end



%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Campbell Diagram (initially x was speed in rpm, now it is in Hz)

figure(1)
plot(N_vector,N_vector,'r','LineWidth',1.5)
hold on
plot(N_vector,lambda_campbell(1:N_campbell,1),'b','LineWidth',1.5)
hold on
plot(N_vector,lambda_campbell(1:N_campbell,3),'b','LineWidth',1.5)
hold on
plot(N_vector,lambda_campbell(1:N_campbell,5),'b','LineWidth',1.5)
hold on
plot(N_vector,lambda_campbell(1:N_campbell,7),'b','LineWidth',1.5)
hold on
plot(N_vector,lambda_campbell(1:N_campbell,9),'b','LineWidth',1.5)
hold on
plot(N_vector,lambda_campbell(1:N_campbell,11),'b','LineWidth',1.5)
hold on
title('Campbell´s Diagram - Critical Speeds for System II','FontSize',16)
xlabel('Angular Velocity N [Hz]','FontSize',14)
ylabel('Imaginary Parts of \lambda_i [Hz]','FontSize',14)
axis([0 450 0 450])
grid on
legend('Unbalance', 'Location', 'northeast')

% Stability Analysis

figure(2)
plot(Omegarpm,lambda_stability(1:N_campbell,1),'b','LineWidth',1.5)
line=yline(0,'--r','LineWidth',1.5);
hold on                                   % this one is the critical one
plot(Omegarpm,lambda_stability(1:N_campbell,3),'b','LineWidth',1.5)
hold on
plot(Omegarpm,lambda_stability(1:N_campbell,5),'b','LineWidth',1.5)
hold on
plot(Omegarpm,lambda_stability(1:N_campbell,7),'b','LineWidth',1.5)
hold on
plot(Omegarpm,lambda_stability(1:N_campbell,9),'b','LineWidth',1.5)
hold on
plot(Omegarpm,lambda_stability(1:N_campbell,11),'b','LineWidth',1.5)
hold on
title('Stability Map - System II','FontSize',14)
xlabel('Angular Velocity \Omega [rpm]','FontSize',16)
ylabel('Real Parts of \lambda_i [rad/s]','FontSize',14)
axis([0 8000 -200 100])
grid on
legend(line,'Stability threshold', 'Location', 'northwest')

findLambda(:,1) = lambda_stability(1:N_campbell,3); % lambda R
findLambda(:,2) = lambda_campbell(1:N_campbell,3);  % lambda I
findLambda % PLotting lambdas to find where it changes sign

