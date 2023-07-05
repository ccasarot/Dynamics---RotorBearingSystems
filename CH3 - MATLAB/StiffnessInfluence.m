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
% Project 2 - Chapter 3                                 %
%                                                       %
%               Copenhagen, Spring semester 2023        %
%                                                       %
%                     Christian Casarotto - s223302     %
%                                                       %
% From only shaft model, study of bearing stiffness     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONTENT
% The following file creates a digital twin of the compressor shaft. In
% this digital twin, the shaft is divivded in elements by nodes, and some
% nodes are assigned to where bearings and disks will fit.
% The process of creating the digital twin is the one describel in chapter
% 1 of the report. 
% 
% As the purpose of this file is to only study the influence of the bearing
% stiffness, no disk is added, and only the shaft is considered. Also, as
% the considerd type of bearing is ball bearing damping is considered null.
% Here bearing stiffness is increased from 0 to 10^9, classic value of ball
% bearing stiffness, to highlight the influence in the shaft mode shapes,
% as descibred in the report.

 clear all
 close all
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DEFINITION OF THE STRUCTURE OF THE MODEL   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 NE=20;         % number of shaft elements
 GL = (NE+1)*4; % number of degree of freedom
 ND=2;          % number of discs
 NM=2;          % number of bearings
 CD1=4;         % node - disc 1
 CD2=10;        % node - disc 2
 CMM1=8;        % node - bearing 1
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

 Rd   = 6/100;                               % external radius of the disc [m]
 Ri   = (5/2)/1000;                          % internal radius of the disc [m]
 espD = 1.1/100 ;                            % disc thickness  [m]
 MasD = pi*Rd^2*espD*RAl - pi*Ri^2*espD*RAl; % disc mass (alunimum) [kg]                                 
 Id   = (1/4*Rd^2+1/12*espD^2)*pi*Rd^2*espD*RAl - (1/4*Ri^2+1/12*espD^2)*pi*Ri^2*espD*RAl;
                                             % transversal mass moment of inertia of the disc [Kgm^2]                                       
 Ip   = 1/2*Rd*Rd*(pi*Rd^2*espD*RAl) - 1/2*Ri*Ri*(pi*Ri^2*espD*RAl);                     
                                             % polar mass moment of inertia of the disc [Kgm^2]
 
 %%%%%%%%%%%%%%%
 %(B) BEARINGS %     % all "Bearings" is later brutally overwritten 
 %%%%%%%%%%%%%%%
 
 MasM = 0.40698;                        % bearing mass [kg](housing + ball bearings)
 h=1/1000;                              % beam thickness [m]
 b=28.5/1000;                           % beam width [m]
 Area=b*h;                              % beam cross section area [m^2] 
 I=b*h^3/12;                            % beam moment of inertia of area [m^4]  
 lr=7.5/100;                            % beam length [m]
 Kty0=2*12*E*I/(lr^3);                  % equivalent beam flexural stiffness [N/m]
 Ktz0=2*E*Area/lr;                      % equivalent bar stiffness [N/m]
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

% length and radius of the shaft elements [m]
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
    for i=1:NE                 ri(i)=Rint; end
% density of shaft elements [kg/m]
    for i=1:NE                 ro(i) = RAco; end
% transversal areal of the shaft elements [m^2]}
    for i=1:NE St(i) = pi*(rx(i)^2-ri(i)^2); end
% area moment of inertia of the shaft elements [m^4]}
    for i=1:NE II(i)=pi*(rx(i)^4-ri(i)^4)/4; end
    


MatrixIndex = 0;

salti_1 = 100;
salti_2 = 1000;
salti_3 = 100000;
salti_4 = 1000000;

for number = [10:salti_1:1000 , 1000:salti_2:10000 , 10000:salti_3:100000 , 100000:salti_4:10^9]

MatrixIndex = MatrixIndex+1;

Stiffness = number; % stiffness of the bearings
StiffnessVector(MatrixIndex) = Stiffness;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       MOUNTING THE GLOBAL MATRICES           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defining the global matrices with zero elements
   M=zeros(GL);
   G=zeros(GL);
   K=zeros(GL);   
   
%%%%%%%%%%%%%%%%%%%%%%   
% GLOBAL MASS MATRIX %
%%%%%%%%%%%%%%%%%%%%%%

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

%    M((CD1-1)*4+1,(CD1-1)*4+1)=M((CD1-1)*4+1,(CD1-1)*4+1)+MasD;
%    M((CD1-1)*4+2,(CD1-1)*4+2)=M((CD1-1)*4+2,(CD1-1)*4+2)+MasD;
%    M((CD1-1)*4+3,(CD1-1)*4+3)=M((CD1-1)*4+3,(CD1-1)*4+3)+Id;
%    M((CD1-1)*4+4,(CD1-1)*4+4)=M((CD1-1)*4+4,(CD1-1)*4+4)+Id;
%    
%    M((CD2-1)*4+1,(CD2-1)*4+1)=M((CD2-1)*4+1,(CD2-1)*4+1)+MasD;
%    M((CD2-1)*4+2,(CD2-1)*4+2)=M((CD2-1)*4+2,(CD2-1)*4+2)+MasD;
%    M((CD2-1)*4+3,(CD2-1)*4+3)=M((CD2-1)*4+3,(CD2-1)*4+3)+Id;
%    M((CD2-1)*4+4,(CD2-1)*4+4)=M((CD2-1)*4+4,(CD2-1)*4+4)+Id;
% 
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

%    G((CD1-1)*4+3,(CD1-1)*4+4)=G((CD1-1)*4+3,(CD1-1)*4+4)-Ip;
%    G((CD1-1)*4+4,(CD1-1)*4+3)=G((CD1-1)*4+4,(CD1-1)*4+3)+Ip;
%    
%    G((CD2-1)*4+3,(CD2-1)*4+4)=G((CD2-1)*4+3,(CD2-1)*4+4)-Ip;
%    G((CD2-1)*4+4,(CD2-1)*4+3)=G((CD2-1)*4+4,(CD2-1)*4+3)+Ip;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL STIFFNESS MATRIX %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
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

 % Part that I added
 NatFreq = diag(lambda);
 NatFreq = NatFreq(1:2:end);
 NatFreq = imag(NatFreq);
 NatFreq = NatFreq/(2*pi);
 NatFreq = sort(NatFreq);

 [lam,p]=sort(diag(lambda));
 U=U(:,p);

 % Number of divisions in time for plotting the mode shapes
  nn=99;

  N=size(U,1);
  maximo=num2str((N-2)/2);
  ModoVirt=N;

  ModoVirt = 1;
  %For visualizing the mode shapes:
  ModoReal=2*ModoVirt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        LOOP TO PLOT THE MODES SHAPES            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 choices = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 16 0];
 orientation = zeros(length(NatFreq),1);
 index = 1;
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

  if all(v==0)
      %fprintf('Vertical \n')
      orientation(index)=1;
  else
      %fprintf('Horizontal \n')
      orientation(index)=0;
  end
  
  Zero=diag(zeros(length(t)))';
  Um=diag(eye(length(t)))';
  nm=num2str(ModoVirt);
  fn=num2str(abs(wn/2/pi));
  dfi=num2str(Omegarpm);
  title(['Angular Velocity (\Omega): ',dfi,' rpm       Mode: ',nm,'   Nat. Freq.:  ',fn,' Hz'],'FontSize',14)
  view(-25,20);
  grid on;
  %ModoVirt=input([' Enter the number of the mode shape to be plotted, zero to esc, highest mode ',maximo,': ']);
  ModoVirt = choices(index);
  index = index + 1;
  if ModoVirt == 0
     break
  end
  ModoReal=2*ModoVirt;
  %figure(ModoVirt+number)
 
 end
scelte = NatFreq(orientation==1);
scelte = scelte(1:8);
MatriceScelte(:,MatrixIndex) = scelte;
end



clear x; clear y;
% Generate x values
NumberOfThicks = size(MatriceScelte, 2);
x = 1:NumberOfThicks;



% to plot values use StiffnessVector
loglog(StiffnessVector, MatriceScelte);
% Set axis labels
xlabel('Bearing Stiffness [N/m]','FontSize',14)
ylabel('Natural frequency [Hz]','FontSize',14)
title('Correlation between natural frequencies and stiffness','FontSize',16)
legend; fontsize(legend,14,'points')



figure
plot(StiffnessVector, MatriceScelte);
% Set axis labels
xlabel('Bearing Stiffness [N/m]','FontSize',14)
ylabel('Natural frequency [Hz]','FontSize',14)
title('Correlation between natural frequencies and stiffness','FontSize',16)


