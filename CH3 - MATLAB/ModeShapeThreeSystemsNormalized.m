%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project 2 - Chapter 3                                 %
%                                                       %
%               Copenhagen, Spring semester 2023        %
%                                                       %
%                     Christian Casarotto - s223302     %
%                                                       %
% Adding rotation, study of critical speeds             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONTENT
% The following file creates a digital twin of the compressor shaft. In
% this digital twin, the shaft is divivded in elements by nodes, and some
% nodes are assigned to where bearings and disks will fit.
% The process of creating the digital twin is the one describel in chapter
% 1 of the report. 
% 
% This file plots the mode shapes corresponding to each system for the
% given bearing stiffness of K=10^9.
% The plot is normalized.

clear all
close all

 Stiffness = 10^9;           % Bearing Stiffness [N/m]

 Omega=  0*2*pi;  % angular velocity [rad/s]
 Omegarpm = Omega*60/2/pi;

 E    = 2.0E11;  % {elasticity modulus [N/m^2}
 RAco = 7800;    % {steel density [kg/m^3]}
 RAl  = 2770;    % {aluminum density [kg/m^3]}



for indice=1:3
indiceSTR = num2str(indice);



if indice==1

 Omega = 263.2; % [Hz]
 Omegarpm = Omega*60;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DEFINITION OF THE STRUCTURE OF THE MODEL   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 NE=20;         % number of shaft elements
 GL = (NE+1)*4; % number of degree of freedom
 ND=2;          % number of discs
 NM=2;          % number of bearings
 CD1=0;         % node - disc 1
 CD2=0;        % node - disc 2
 CMM1=8;        % node - bearing 1
 CMM2=NE;       % node - bearing 2

 % (C) SHAFT  
 Rint = 0;                     % shaft internal radius [m]                      

 l(1)  = 65/1000;    rx(1)  = 54.28/2000;
 l(2)  = 35/1000;    rx(2)  = 70/2000;
 l(3)  = 40/1000;    rx(3)  = 70/2000;
 l(4)  = 50/1000;    rx(4)  = 70/2000;
 l(5)  = 50/1000;    rx(5)  = 70/2000;
 l(6)  = 55/1000;    rx(6)  = 90/2000;
 l(7)  = 60/1000;    rx(7)  = 99.6/2000;
 l(8)  = 60/1000;    rx(8)  = 99.6/2000;
 l(9)  = 62.5/1000;  rx(9)  = 90/2000;
 l(10) = 62.5/1000;  rx(10) = 90/2000;
 l(11) = 62.5/1000;  rx(11) = 90/2000;
 l(12) = 62.5/1000;  rx(12) = 90/2000;
 l(13) = 62.5/1000;  rx(13) = 90/2000;
 l(14) = 62.5/1000;  rx(14) = 90/2000;
 l(15) = 62.5/1000;  rx(15) = 90/2000;
 l(16) = 62.5/1000;  rx(16) = 90/2000;
 l(17) = 62.5/1000;  rx(17) = 90/2000;
 l(18) = 62.5/1000;  rx(18) = 90/2000;
 l(19) = 46.5/1000;  rx(19) = 53.27/2000; 
 l(20) = 63.5/1000;  rx(20) = 53.27/2000; 

% internal radius of shaft elements [m]
    for i=1:NE                        ri(i)=Rint; end
% density of shaft elements [kg/m]
    for i=1:NE                      ro(i) = RAco; end
% transversal areal of the shaft elements [m^2]}
    for i=1:NE      St(i) = pi*(rx(i)^2-ri(i)^2); end
% area moment of inertia of the shaft elements [m^4]}
    for i=1:NE      II(i)=pi*(rx(i)^4-ri(i)^4)/4; end
    
% MOUNTING THE GLOBAL MATRICES     

   M=zeros(GL); G=zeros(GL); K=zeros(GL);   
     
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

% Adding the stiffness matrices of the bearing elements

   K((CMM1-1)*4+1,(CMM1-1)*4+1)=K((CMM1-1)*4+1,(CMM1-1)*4+1)+ Stiffness;
   K((CMM1-1)*4+2,(CMM1-1)*4+2)=K((CMM1-1)*4+2,(CMM1-1)*4+2)+ Stiffness;
   K((CMM2-1)*4+1,(CMM2-1)*4+1)=K((CMM2-1)*4+1,(CMM2-1)*4+1)+ Stiffness;
   K((CMM2-1)*4+2,(CMM2-1)*4+2)=K((CMM2-1)*4+2,(CMM2-1)*4+2)+ Stiffness;

end

if indice==2

 Omega = 133.5; % [Hz]
 Omegarpm = Omega*60;

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
 CMM2=NE;       % node - bearing 2
 
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
 Rint = 0;       % shaft internal radius [m]                      

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
    for i=1:NE                   ri(i)=Rint; end
% density of shaft elements [kg/m]
    for i=1:NE                 ro(i) = RAco; end
% transversal areal of the shaft elements [m^2]}
    for i=1:NE St(i) = pi*(rx(i)^2-ri(i)^2); end
% area moment of inertia of the shaft elements [m^4]}
    for i=1:NE II(i)=pi*(rx(i)^4-ri(i)^4)/4; end
    
% MOUNTING THE GLOBAL MATRICES          
M=zeros(GL); G=zeros(GL); K=zeros(GL);   
   
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
   M((CD2-1)*4+3,(CD2-1)*4+3)=M((CD2-1)*4+3,(CD2-1)*4+3)+Id_2;
   M((CD2-1)*4+4,(CD2-1)*4+4)=M((CD2-1)*4+4,(CD2-1)*4+4)+Id_2;
               
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

% Adding the stiffness matrices of the bearing elements

   K((CMM1-1)*4+1,(CMM1-1)*4+1)=K((CMM1-1)*4+1,(CMM1-1)*4+1) + Stiffness;
   K((CMM1-1)*4+2,(CMM1-1)*4+2)=K((CMM1-1)*4+2,(CMM1-1)*4+2) + Stiffness;
   K((CMM2-1)*4+1,(CMM2-1)*4+1)=K((CMM2-1)*4+1,(CMM2-1)*4+1) + Stiffness;
   K((CMM2-1)*4+2,(CMM2-1)*4+2)=K((CMM2-1)*4+2,(CMM2-1)*4+2) + Stiffness;

end

if indice==3

 Omega = 89.1; % [Hz]
 Omegarpm = Omega*60;

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
 CMM2=NE;       % node - bearing 2
 
% (A) DISCS   

% ------------------ DISK 1
 Rd_1   = 295/2000;          % external radius of the disc [m]
 Ri_1   = 70/2000;        % internal radius of the disc [m]
 espD_1 = 80/1000 ;         % disc thickness  [m]
 % Removed ring
 Rr_1 = 271/2000;             % Internal ratius ring
 Lr_1 = (35/3)/1000;          % Lenght of the ring
 Rdist_1 = 25/1000;           % Ring distance
 % Other rings
 R_external_SmallRing_1 = 90/2000;             % Internal ratius ring
 L_SmallRing_1 = 20/1000;          % Lenght of the ring
 Rdist_SmallRing_1 = 5/1000;           % Ring distance
 R_external_BigRing_1 = 110/2000;             % Internal ratius ring
 L_BigRing_1 = 25/1000;          % Lenght of the ring
 Rdist_BigRing_1 = 27.5/1000;           % Ring distance

 % disc mass [kg]  
 Mde_1 = pi*Rd_1^2*espD_1*RAco;  % Mass disc external
 Mdi_1 = pi*Ri_1^2*espD_1*RAco;  % Mass disc internal
 Mre_1 = pi*Rd_1^2*Lr_1*RAco;    % Mass ring external
 Mri_1 = pi*Rr_1^2*Lr_1*RAco;    % Mass ring internal
     M_external_SmallRing_1 = pi*R_external_SmallRing_1^2*L_SmallRing_1*RAco;    
     M_internal_SmallRing_1 = pi*Ri_1^2*L_SmallRing_1*RAco;    
     M_external_BigRing_1 = pi*R_external_BigRing_1^2*L_BigRing_1*RAco;    
     M_internal_BigRing_1 = pi*Ri_1^2*L_BigRing_1*RAco;     
 MasD_1 = Mde_1 - Mdi_1 - 3*(Mre_1 - Mri_1) - ...
     (M_external_SmallRing_1 - M_internal_SmallRing_1) - ...
     (M_external_BigRing_1 - M_internal_BigRing_1);

 % Inertia [Kgm^2]
 Id_1 = ((1/4)*Rd_1^2 + (1/12)*espD_1^2) * Mde_1 ... % Overall ext
       -((1/4)*Ri_1^2 + (1/12)*espD_1^2) * Mdi_1 ... % Overall int
       -((1/4)*Rd_1^2 + (1/12)*Lr_1^2) * Mre_1 ...   % Central ring ext
       +((1/4)*Rr_1^2 + (1/12)*Lr_1^2) * Mri_1 ...   % Central ring int
       - 2 * ((1/4)*Rd_1^2 + (1/12)*Lr_1^2) * Mre_1 ...   % Side rings ext
       + 2 * ((1/4)*Rr_1^2 + (1/12)*Lr_1^2) * Mri_1 ...   % Side rings int
       - 2 * (Mre_1-Mri_1) * Rdist_1^2 ...                   % Transport 
       - ((1/4)*R_external_SmallRing_1^2 + (1/12)*L_SmallRing_1^2) * M_external_SmallRing_1 ...   
       + ((1/4)*Ri_1^2 + (1/12)*L_SmallRing_1^2) * M_internal_SmallRing_1 ...   
       - (M_external_SmallRing_1-M_internal_SmallRing_1) * Rdist_SmallRing_1^2 ...
       - ((1/4)*R_external_BigRing_1^2 + (1/12)*L_BigRing_1^2) * M_external_BigRing_1 ...   
       + ((1/4)*Ri_1^2 + (1/12)*L_BigRing_1^2) * M_internal_BigRing_1 ...   
       - (M_external_BigRing_1-M_internal_BigRing_1) * Rdist_BigRing_1^2;  

 Ip_1 =   1/2*(Rd_1^2)*Mde_1 ... % Overall disk
        - 1/2*(Ri_1^2)*Mdi_1 ...
            - 1/2*(Rd_1^2)*Mre_1 ... % External ring (scanalature)
            + 1/2*(Rr_1^2)*Mri_1 ...
        - 1/2*(R_external_SmallRing_1^2)*M_external_SmallRing_1 ... % Internal ring piccolo (sedi)
        + 1/2*(Ri_1^2)*M_internal_SmallRing_1 ...
            - 1/2*(R_external_BigRing_1^2)*M_external_BigRing_1 ... % Internal ring grande (sedi)
            + 1/2*(Ri_1^2)*M_internal_BigRing_1;

% ------------------ DISK 2
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
 Rint = 0;       % shaft internal radius [m]                      

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
    for i=1:NE                   ri(i)=Rint; end
% density of shaft elements [kg/m]
    for i=1:NE                 ro(i) = RAco; end
% transversal areal of the shaft elements [m^2]}
    for i=1:NE St(i) = pi*(rx(i)^2-ri(i)^2); end
% area moment of inertia of the shaft elements [m^4]}
    for i=1:NE II(i)=pi*(rx(i)^4-ri(i)^4)/4; end
    
% MOUNTING THE GLOBAL MATRICES          
M=zeros(GL); G=zeros(GL); K=zeros(GL);     

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

% Adding the mass matrices of the disc elements 

   M((CD1-1)*4+1,(CD1-1)*4+1)=M((CD1-1)*4+1,(CD1-1)*4+1)+MasD_1;
   M((CD1-1)*4+2,(CD1-1)*4+2)=M((CD1-1)*4+2,(CD1-1)*4+2)+MasD_1;
   M((CD1-1)*4+3,(CD1-1)*4+3)=M((CD1-1)*4+3,(CD1-1)*4+3)+Id_1*1.1;
   M((CD1-1)*4+4,(CD1-1)*4+4)=M((CD1-1)*4+4,(CD1-1)*4+4)+Id_1*1.1;
   
   M((CD2-1)*4+1,(CD2-1)*4+1)=M((CD2-1)*4+1,(CD2-1)*4+1)+MasD_2;
   M((CD2-1)*4+2,(CD2-1)*4+2)=M((CD2-1)*4+2,(CD2-1)*4+2)+MasD_2;
   M((CD2-1)*4+3,(CD2-1)*4+3)=M((CD2-1)*4+3,(CD2-1)*4+3)+Id_2*1.1;
   M((CD2-1)*4+4,(CD2-1)*4+4)=M((CD2-1)*4+4,(CD2-1)*4+4)+Id_2*1.1;
           
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

% Adding the gyroscopic matrices of the disc elements

   G((CD1-1)*4+3,(CD1-1)*4+4)=G((CD1-1)*4+3,(CD1-1)*4+4)-Ip_1;
   G((CD1-1)*4+4,(CD1-1)*4+3)=G((CD1-1)*4+4,(CD1-1)*4+3)+Ip_1;
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

% Adding the stiffness matrices of the bearing elements
   K((CMM1-1)*4+1,(CMM1-1)*4+1)=K((CMM1-1)*4+1,(CMM1-1)*4+1) + Stiffness;
   K((CMM1-1)*4+2,(CMM1-1)*4+2)=K((CMM1-1)*4+2,(CMM1-1)*4+2) + Stiffness;
   K((CMM2-1)*4+1,(CMM2-1)*4+1)=K((CMM2-1)*4+1,(CMM2-1)*4+1) + Stiffness;
   K((CMM2-1)*4+2,(CMM2-1)*4+2)=K((CMM2-1)*4+2,(CMM2-1)*4+2) + Stiffness;


end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %    GLOBAL MATHEMATICAL MODEL                 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
  Mglob=[   M                zeros(size(M,1)) 
            zeros(size(M,1)) K             ];
  Kglob=[  -Omega*G          K 
           -K                zeros(size(M,1))];
      
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %              MODAL ANALYSIS                  %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 % Calculating Eigenvectors and Eigenvalues

 [U,lambda]=eig(-Kglob,Mglob);
 [lam,p]=sort(diag(lambda));
 U=U(:,p);

 % Number of divisions in time for plotting the mode shapes
  nn=99;

  N=size(U,1);
  maximo=num2str((N-2)/2);
  ModoVirt=N;

  ModoVirt=input(['Step 1 indice ',indiceSTR,': ']);

  %For visualizing the mode shapes:

  ModoReal=2*ModoVirt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        LOOP TO PLOT THE MODES SHAPES            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 while ModoReal>0

  % Natural frequencies
  wn=imag(lam(ModoReal));
  ttotal=8/abs(wn);
  dt=ttotal/nn;
  t=0:dt:ttotal;

  % Defining v and w real e imaginary for each node

  y=1:4:GL;
  z=2:4:GL;

  for i=1:(NE+1)
   vr(i)=real(U(y(i),ModoReal));
   vi(i)=imag(U(y(i),ModoReal));
   wr(i)=real(U(z(i),ModoReal));
   wi(i)=imag(U(z(i),ModoReal));
  end

  % Calculation the modal displacements v and w

  for i=1:(NE+1)
   v(i,:)=vr(i)*cos(wn*t)+vi(i)*sin(wn*t);
   w(i,:)=wr(i)*cos(wn*t)+wi(i)*sin(wn*t);
  end
  
  Zero=diag(zeros(length(t)))';
  Um=diag(eye(length(t)))';

  for i=1:(NE+1)
   pos(i,:)=Zero+(i-1)*Um;
  end
 

if indice ==1 
  % axis scale for plotting the mode shapes  
  scale_VW = max(max(abs(v))) + max(max(abs(w)));
  scale_pos = max(max(abs(pos)));
end

  %clf
  hold on
  for cont=1:NE+1
        if indice==1
            plot_v_1 = v;
            plot_w_1 = w;
            plot_v = plot_v_1;
            plot_w = plot_w_1;
        end
        if indice==2
            plot_v_2 = v;
            plot_w_2 = w;
            plot_v = plot_v_2;
            plot_w = plot_w_2;
        end
        if indice==3
            plot_v_3 = v;
            plot_w_3 = w;
            plot_v = plot_v_3;
            plot_w = plot_w_3;
        end
  end

  ModoVirt=input(['Step 2 indice ',indiceSTR,': ']);
  if ModoVirt==50
      break
  end
  ModoReal=2*ModoVirt;
 
 end

 % END MAIN LOOP
end



indice=0;
for j=1:3
    indice=indice+1;

if indice ==1 color = [230, 159, 0] / 255; 
    correction=0;
            plot_v = plot_v_1;
            plot_w = plot_w_1;
            max1=max(plot_v(13,:)); % if you do w its the same
            plot_v=plot_v/max1;
            plot_w=plot_w/max1;
end % red
if indice ==2 color = [0, 158, 115] / 255; 
    correction=0.2;
            plot_v = plot_v_2;
            plot_w = plot_w_2;
            max2=max(plot_v(13,:)); % if you do w its the same
            plot_v=plot_v/max2;
            plot_w=plot_w/max2;
end % green
if indice ==3 color = [0, 114, 178] / 255; 
    correction=-0.2;
            plot_v = plot_v_3;
            plot_w = plot_w_3;
            max3=max(plot_v(13,:)); % if you do w its the same
            plot_v=plot_v/max3;
            plot_w=plot_w/max3;
end

NE = 20;
for cont=1:NE+1
    hold on
    plot3(pos(cont,:)+correction,plot_w(cont,:),plot_v(cont,:),'Color',[color,0.7],'LineWidth',3);
    axis([-0.2,scale_pos+0.2,-3,3,-3,3])
    title('Mode shape 1 for the three systems','FontSize',16)
    xlabel('Shaft nodes','FontSize',14)
    ylabel('Modal displ. "v"','FontSize',14)
    zlabel('Modal displ. "w"','FontSize',14)
    view(-25,20);
    grid on;
end
end