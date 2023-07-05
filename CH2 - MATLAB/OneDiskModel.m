% Original file header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MACHINERY DYNAMICS LECTURES  (41514)                  %
% MEK - DEPARTMENT OF MECHANICAL ENGINEERING            %
% DTU - TECHNICAL UNIVERSITY OF DENMARK                 %
%                                                       %
%              Copenhagen, March 30th, 2021             %
%                                                       %
%                        Ilmar Ferreira Santos          %
%                                                       %
% ROTATING MACHINES -- NATURAL FREQUENCIES AND MODES    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Modified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project 2 - Chapter 2                                 %
%                                                       %
%               Copenhagen, Spring semester 2023        %
%                                                       %
%                     Christian Casarotto - s223302     %
%                                                       %
% Adding one disk to the improved model                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONTENT
% The following file creates a digital twin of the compressor shaft. In
% this digital twin, the shaft is divivded in elements by nodes, and some
% nodes are assigned to where bearings and disks will fit.
% The process of creating the digital twin is the one describel in chapter
% 1 of the report.
% The bearings are of course the supports of the shaft. Still no damping is
% given to the bearings, beside a little amount to avoid numerical
% problems.
% The file runs a loop where the user needs to select the natural frequancy
% of which he wants to plot mode shape and value.
% 
% Here the bigger of the two compressor disks is added to the model.
% Several versions of the file were created, with different grades of
% precision in the modelling of the disk based on the technical drawings.
% Finally the best solution, giving results close to the experimental data,
% was found for a very detailed discretization of the disk, as described in
% the report. Also, a reduction factor of the disk inertia was added to
% better fit the data.

 clear all
 close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DEFINITION OF THE STRUCTURE OF THE MODEL   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 NE=22;         % number of shaft elements
 GL = (NE+1)*4; % number of degree of freedom
 ND=2;          % number of discs
 NM=2;          % number of bearings
 CD1=4;         % node - disc 1
 CD2=6;         % node - disc 2
 CMM1=1;        % node - bearing 1
 CMM2=NE;       % node - bearing 2
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               CONSTANTS                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 E    = 2.0E11;  % {elasticity modulus [N/m^2}
 RAco = 7800;    % {steel density [kg/m^3]}
 RAl  = 2770;    % {aluminum density [kg/m^3]}
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         OPERATIONAL CONDITIONS               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Omega=  0*2*pi;  % angular velocity [rad/s]
  Omegarpm = Omega*60/2/pi; 
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   GEOMETRY OF THE ROTATING MACHINE           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%
%(A) DISCS    %
%%%%%%%%%%%%%%%

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

 % Inertia
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

                             % transversal mass moment of inertia [Kgm^2]                                       
 Ip_2   = 1/2*(Rd_2^2)*(pi*Rd_2^2*espD_2*RAco) ...
     - 1/2*(Ri_2^2)*(pi*Ri_2^2*espD_2*RAco) ...
     + 1/2*(Rd_2^2)*(pi*Rd_2^2*Lr_2*RAco) ...
     - 1/2*(Rr_2^2)*(pi*Rr_2^2*Lr_2*RAco); % - pi*((Rd_2)^2*Lr_2-(Rr_2)^2*Lr_2)*RAco*Rr_2^2                  
                             % polar mass moment of inertia [Kgm^2]
 
 %%%%%%%%%%%%%%%
 %(B) BEARINGS %     % all "Bearings" is later brutally overwritten 
 %%%%%%%%%%%%%%%
 
 MasM = 0.40698;                      % bearing mass [kg](housing + ball bearings)
 h=1/1000;                            % beam thickness [m]
 b=28.5/1000;                         % beam width [m]
 Area=b*h;                            % beam cross section area [m^2] 
 I=b*h^3/12;                          % beam moment of inertia of area [m^4]  
 lr=7.5/100;                          % beam length [m]
 Kty0=2*12*E*I/(lr^3);                % equivalent beam flexural stiffness [N/m]
 Ktz0=2*E*Area/lr;                    % equivalent bar stiffness [N/m]
% Bearing 1 - Damping
 Dty1 = 0.0 ;                 
 Dtz1 = 0.0 ;                 
 Dry1 = 0.0 ;                
 Drz1 = 0.0 ;                  
% Bearing 2 - Damping
 Dty2 = 0.0 ;                 
 Dtz2 = 0.0 ;                 
 Dry2 = 0.0 ;                 
 Drz2 = 0.0 ;  
% Bearing 1 - Stiffness
 Kty1 = Kty0 ;        
 Ktz1 = Ktz0 ;            
 Kry1 = 0.0 ;                 
 Krz1 = 0.0 ;                 
% Bearing 2 - Stiffness
 Kty2 = Kty0 ;                
 Ktz2 = Ktz0 ;                
 Kry2 = 0.0 ;                 
 Krz2 = 0.0 ;  
 
 %%%%%%%%%%%%%%%
 %(C) SHAFT    %
 %%%%%%%%%%%%%%%
 
 Rint = 0;       % shaft internal radius [m]                      

 % length of the shaft elements [m]
 % external radius of shaft elements [m]

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
 l(21) = 23/1000;      rx(21) = 50/2000; 
 l(22) = 52/1000;      rx(22) = 40/2000; 

% internal radius of shaft elements [m]
    for i=1:NE                   ri(i)=Rint; end
% density of shaft elements [kg/m]
    for i=1:NE                 ro(i) = RAco; end
% transversal areal of the shaft elements [m^2]}
    for i=1:NE St(i) = pi*(rx(i)^2-ri(i)^2); end
% area moment of inertia of the shaft elements [m^4]}
    for i=1:NE II(i)=pi*(rx(i)^4-ri(i)^4)/4; end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       MOUNTING THE GLOBAL MATRICES           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 disp('MOUNTING THE GLOBAL MATRICES - WAIT!')
 disp(' ')

% Defining the global matrices with zero elements
   M=zeros(GL);
   G=zeros(GL);
   K=zeros(GL);   
   
%%%%%%%%%%%%%%%%%%%%%%   
% GLOBAL MASS MATRIX %
%%%%%%%%%%%%%%%%%%%%%%

disp('MOUNTING THE GLOBAL MASS MATRIX - WAIT!')
disp(' ')

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

%    M((CD1-1)*4+1,(CD1-1)*4+1)=M((CD1-1)*4+1,(CD1-1)*4+1)+MasD_1;
%    M((CD1-1)*4+2,(CD1-1)*4+2)=M((CD1-1)*4+2,(CD1-1)*4+2)+MasD_1;
%    M((CD1-1)*4+3,(CD1-1)*4+3)=M((CD1-1)*4+3,(CD1-1)*4+3)+Id_1;
%    M((CD1-1)*4+4,(CD1-1)*4+4)=M((CD1-1)*4+4,(CD1-1)*4+4)+Id_1;
   
InertiaFactor = 0.9;

   M((CD2-1)*4+1,(CD2-1)*4+1)=M((CD2-1)*4+1,(CD2-1)*4+1)+MasD_2;
   M((CD2-1)*4+2,(CD2-1)*4+2)=M((CD2-1)*4+2,(CD2-1)*4+2)+MasD_2;
   M((CD2-1)*4+3,(CD2-1)*4+3)=M((CD2-1)*4+3,(CD2-1)*4+3)+Id_2*InertiaFactor;
   M((CD2-1)*4+4,(CD2-1)*4+4)=M((CD2-1)*4+4,(CD2-1)*4+4)+Id_2*InertiaFactor;

% % Adding the mass matrices of the bearing elements 
% 
%    M((CMM1-1)*4+1,(CMM1-1)*4+1)=M((CMM1-1)*4+1,(CMM1-1)*4+1)+MasM;
%    M((CMM1-1)*4+2,(CMM1-1)*4+2)=M((CMM1-1)*4+2,(CMM1-1)*4+2)+MasM;
%    
%    M((CMM2-1)*4+1,(CMM2-1)*4+1)=M((CMM2-1)*4+1,(CMM2-1)*4+1)+MasM;
%    M((CMM2-1)*4+2,(CMM2-1)*4+2)=M((CMM2-1)*4+2,(CMM2-1)*4+2)+MasM;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL GYROSCOPIC MATRIX %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('MOUNTING THE GLOBAL GYROSCOPIC MATRIX - WAIT!')
disp(' ')
           
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

%    G((CD1-1)*4+3,(CD1-1)*4+4)=G((CD1-1)*4+3,(CD1-1)*4+4)-Ip_2;
%    G((CD1-1)*4+4,(CD1-1)*4+3)=G((CD1-1)*4+4,(CD1-1)*4+3)+Ip_2;
   
   G((CD2-1)*4+3,(CD2-1)*4+4)=G((CD2-1)*4+3,(CD2-1)*4+4)-Ip_2;
   G((CD2-1)*4+4,(CD2-1)*4+3)=G((CD2-1)*4+4,(CD2-1)*4+3)+Ip_2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL STIFFNESS MATRIX %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('MOUNTING THE GLOBAL STIFFNESS MATRIX - WAIT!')
disp(' ')
           
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

% REDUCTION FACTOR OF THE BEARING STIFFNESS to avoid numerical problems
RedFactor = 0.0000001;


% Adding the stiffness matrices of the bearing elements

   K((CMM1-1)*4+1,(CMM1-1)*4+1)=K((CMM1-1)*4+1,(CMM1-1)*4+1)+RedFactor*Ktz1;
   K((CMM1-1)*4+2,(CMM1-1)*4+2)=K((CMM1-1)*4+2,(CMM1-1)*4+2)+RedFactor*Kty1;
   
   K((CMM2-1)*4+1,(CMM2-1)*4+1)=K((CMM2-1)*4+1,(CMM2-1)*4+1)+RedFactor*Ktz2;
   K((CMM2-1)*4+2,(CMM2-1)*4+2)=K((CMM2-1)*4+2,(CMM2-1)*4+2)+RedFactor*Kty2;

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
 
 disp('CALCULATING NATURAL FREQUENCIES AND MODE SHAPES - WAIT!')
 disp(' ')

 % Calculating Eigenvectors and Eigenvalues

 [U,lambda]=eig(-Kglob,Mglob);
 [lam,p]=sort(diag(lambda));
 U=U(:,p);

 % Number of divisions in time for plotting the mode shapes
  nn=99;

  N=size(U,1);
  maximo=num2str((N-2)/2);
  ModoVirt=N;

  ModoVirt=input([' Enter the number of the mode shape to be plotted, zero to esc, highest mode ',maximo,': ']);

  %For visualizing the mode shapes:

  ModoReal=2*ModoVirt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        LOOP TO PLOT THE MODES SHAPES        %
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
  
  % axis scale for plotting the mode shapes  
  scale_VW = max(max(abs(v))) + max(max(abs(w)))
  scale_pos = max(max(abs(pos)))
  
  clf
  hold on
  for cont=1:NE+1
   plot3(pos(cont,:),w(cont,:),v(cont,:),'k','LineWidth',2.5);
   axis([0,scale_pos,-scale_VW,scale_VW,-scale_VW, scale_VW])
  end

  nm=num2str(ModoVirt);
  fn=num2str(abs(wn/2/pi));
  dfi=num2str(Omegarpm);
  title(['Angular Velocity: ',dfi,' rpm       Mode: ',nm,'   Nat. Freq.:  ',fn,' Hz'],'FontSize',14)
  view(-25,20);
  grid on;
  ModoVirt=input([' Enter the number of the mode shape to be plotted, zero to esc, highest mode ',maximo,': ']);
  ModoReal=2*ModoVirt;
  figure(ModoVirt)
 
end