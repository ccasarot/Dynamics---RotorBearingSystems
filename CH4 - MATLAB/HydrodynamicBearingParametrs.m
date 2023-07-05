%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project 2 - Chapter 4                                 %
%                                                       %
%               Copenhagen, Spring semester 2023        %
%                                                       %
%                     Christian Casarotto - s223302     %
%                                                       %
% Hydrodynamic bearing parameters                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CONTENT
% In this file the bearing parameters are calculated and plotted as
% described in the report

close all
clear all



% Definition of Nondimensional Bearing Parameters
W = 892.7;                              % W - external load in N; 
eta = 0.0277 * exp(1)^(0.034*(40-55)); % viscosity % η - oil viscosity

% Cross section
% d and r for the shaft D and R for the bearing recall C = R - r

C = 1e-4;                              % [m] clearance % C is the bearing clearance in m; 
d = 99.6/1000;                         % [m] diameter of the shaft % D is the bearing inner diameter in m;
r = d/2; 
R = C + r; % form C = R - r
D = R*2;   % from R = D/2;
syms L
L = solve(L/D==0.5,L);                 % as L/d is 0.5 % L is the bearing width in m; 

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
 
% Journal Bearing -- Static Properties   
% figure(1)
% subplot(2,2,1), plot(Table(:,1),Table(:,2),'*-b','LineWidth',1.5)
% title('Excentricity','FontSize',14)
% xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
% ylabel('\epsilon=e/C','FontSize',14)
% grid    
% subplot(2,2,2), plot(Table(:,1),Table(:,3),'*-b','LineWidth',1.5)
% title('Atittude Angle','FontSize',14)
% xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
% ylabel('\phi [^o]','FontSize',14)
% grid
% subplot(2,2,3.5),

% Locus plot
polar(3*pi/2+Table(:,3)*pi/180,Table(:,2))
title('Journal Center Locus','FontSize',14)
grid

% 
% % Journal Bearing -- Dynamic Properties (Stiffness)
% figure(2)
% subplot(2,2,1), plot(Table(:,1),Table(:,7),'*-b','LineWidth',1.5)
% xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
% ylabel('Kxx','FontSize',14)
% grid    
% subplot(2,2,2), plot(Table(:,1),Table(:,8),'*-b','LineWidth',1.5)
% xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
% ylabel('Kxy')
% grid
% subplot(2,2,3), plot(Table(:,1),Table(:,9),'*-b','LineWidth',1.5)
% xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
% ylabel('Kyx','FontSize',14)
% grid
% subplot(2,2,4), plot(Table(:,1),Table(:,10),'*-b','LineWidth',1.5)
% xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
% ylabel('Kyy','FontSize',14)
% grid
% 
% % Journal Bearing -- Dynamic Properties (Damping)
% figure(3)
% subplot(2,2,1), plot(Table(:,1),Table(:,11),'*-b','LineWidth',1.5)
% xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
% ylabel('Bxx','FontSize',14)
% grid
% subplot(2,2,2), plot(Table(:,1),Table(:,12),'*-b','LineWidth',1.5)
% xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
% ylabel('Bxy','FontSize',14)
% grid
% subplot(2,2,3), plot(Table(:,1),Table(:,13),'*-b','LineWidth',1.5)
% xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
% ylabel('Byx','FontSize',14)
% grid
% subplot(2,2,4), plot(Table(:,1),Table(:,14),'*-b','LineWidth',1.5)
% xlabel('S=\eta*N*L*D/W*(R/C)^2','FontSize',14)
% ylabel('Byy','FontSize',14)
% grid



% Building stiffness and damping
x_st = []; y_st = []; LayerThickness = []; N_vector = [];
kxx = []; kxy = []; kyx = []; kyy = [];
dxx = []; dxy = []; dyx = []; dyy = [];
for line=1:length(Table)

    syms N
    S = eta*N*L*D/W*(R/C)^2;                          % S - Sommerfeld Number
    N = double(solve(S == Table(line,1), N));   % find N for given S
    N_vector(line) = N;
    w = 2*pi*N;                                 % ω rotor angular velocity [rad/s]
        
    % Coefficients of the Stiffness Matrix  
    kxx(line) =  Table(line,7)*W/C  ; % [N/m]  
    kxy(line) =  Table(line,8)*W/C  ; % [N/m]
    kyx(line) =  Table(line,9)*W/C  ; % [N/m]
    kyy(line) =  Table(line,10)*W/C ; % [N/m]
    
    % Coefficients of the Damping Matrix 
    dxx(line) =  Table(line,11)*W/(w*C)  ; % [N/(m/s)]
    dxy(line) =  Table(line,12)*W/(w*C)  ; % [N/(m/s)]
    dyx(line) =  Table(line,13)*W/(w*C)  ; % [N/(m/s)]
    dyy(line) =  Table(line,14)*W/(w*C)  ; % [N/(m/s)]
    
    x_st(line) = Table(line,2)*cos(Table(line,3)*pi/180)*C; % equilibrium position [m]
    y_st(line) = Table(line,2)*sin(Table(line,3)*pi/180)*C; % equilibrium position [m]
    LayerThickness(line) = (D/2) - (d/2) - sqrt(x_st(line)^2+y_st(line)^2);
end
LayerThickness=LayerThickness*1000*1000; % convert to micron



% % Plotting GAP
% figure;
% plot(N_vector, LayerThickness, '*-r', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'r');
% title('Film Thickness as a Function of Speed N','FontSize',16)
% xlabel('Velocity N [Hz]','FontSize',14)
% ylabel('Film Thickness [\mum]','FontSize',14)
% % Adding asterisks and connecting lines
% hold on;
% plot(N_vector, LayerThickness, 'k--');
% scatter(N_vector, LayerThickness, 60, 'r', 'filled');
% text(N_vector, LayerThickness, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
% % Adjusting axis limits
% xlim([min(N_vector) - 1, max(N_vector) + 1]);
% ylim([0, 100]);



% Plotting the original data - Film Thickness as a Function of Speed N
figure;
plot(N_vector, LayerThickness, '*-b', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'b');
title('Film Thickness as a Function of Speed N', 'FontSize', 16)
xlabel('Velocity N [Hz]', 'FontSize', 14)
ylabel('Film Thickness [\mum]', 'FontSize', 14)
hold on;
plot(N_vector, LayerThickness, 'r');
text(N_vector, LayerThickness, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
xlim([min(N_vector) - 1, max(N_vector) + 1]);
ylim([0, 100]);

% Interpolation
N_interpol = min(N_vector):0.1:max(N_vector); % New x-values for interpolation
LayerThickness_interpol = interp1(N_vector, LayerThickness, N_interpol, 'spline'); % Interpolation using spline

% Plotting the interpolated data
plot(N_interpol, LayerThickness_interpol, 'r', 'LineWidth', 2);
legend('Original Data', 'Interpolated Data');



% % Plotting K
% figure;
% subplot(2, 2, 1);
% plot(N_vector, kxx, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
% title('kxx [N/m]', 'FontSize', 14);
% xlabel('N [Hz]','FontSize',14)
% ylabel('kxx [N/m]','FontSize',14)
% subplot(2, 2, 2);
% plot(N_vector, kxy, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
% title('kxy [N/m]', 'FontSize', 14);
% xlabel('N [Hz]','FontSize',14)
% ylabel('kxy [N/m]','FontSize',14)
% subplot(2, 2, 3);
% plot(N_vector, kyx, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
% title('kyx [N/m]', 'FontSize', 14);
% xlabel('N [Hz]','FontSize',14)
% ylabel('kyx [N/m]','FontSize',14)
% subplot(2, 2, 4);
% plot(N_vector, kyy, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
% title('kyy [N/m]', 'FontSize', 14);
% xlabel('N [Hz]','FontSize',14)
% ylabel('kyy [N/m]','FontSize',14)
% 
% % Plotting D
% figure;
% subplot(2, 2, 1);
% plot(N_vector, dxx, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
% title('dxx [N/m]', 'FontSize', 14);
% xlabel('N [Hz]','FontSize',14)
% ylabel('dxx [N/m]','FontSize',14)
% subplot(2, 2, 2);
% plot(N_vector, dxy, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
% title('dxy [N/m]', 'FontSize', 14);
% xlabel('N [Hz]','FontSize',14)
% ylabel('dxy [N/m]','FontSize',14)
% subplot(2, 2, 3);
% plot(N_vector, dyx, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
% title('dyx [N/m]', 'FontSize', 14);
% xlabel('N [Hz]','FontSize',14)
% ylabel('dyx [N/m]','FontSize',14)
% subplot(2, 2, 4);
% plot(N_vector, dyy, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
% title('dyy [N/m]', 'FontSize', 14);
% xlabel('N [Hz]','FontSize',14)
% ylabel('dyy [N/m]','FontSize',14)

% % Plotting K - 2
% spacing=5;              % you always get one less
% N_vector_new = [];
% for i=1:length(N_vector)-1
%     gap = (N_vector(i+1) - N_vector(i))/spacing;
%     vector = [];
%     for k=1:spacing
%         vector(k) = N_vector(i) + gap*k;
%     end
%     N_vector_new = [N_vector_new, vector];
% end
% 
% N_vector_new = [N_vector(1), N_vector_new]
% 
% kxx_new = interp1(kxx, N_vector, N_vector_new, 'spline');
% kxy_new = interp1(kxy, N_vector, N_vector_new, 'spline');
% kyx_new = interp1(kyx, N_vector, N_vector_new, 'spline');
% for i=1:length(kyy) % fix matlab madness
%     factor = rand/100000 + 1
%     kyy(i) = kyy(i)*factor
% end
% kyy_new = interp1(kyy, N_vector, N_vector_new, 'spline');
% 
% figure;
% subplot(2, 2, 1);
% plot(N_vector_new, kxx_new, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
% title('kxx [N/m]', 'FontSize', 14);
% xlabel('N [Hz]','FontSize',14)
% ylabel('kxx [N/m]','FontSize',14)
% subplot(2, 2, 2);
% plot(N_vector_new, kxy_new, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
% title('kxy [N/m]', 'FontSize', 14);
% xlabel('N [Hz]','FontSize',14)
% ylabel('kxy [N/m]','FontSize',14)
% subplot(2, 2, 3);
% plot(N_vector_new, kyx_new, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
% title('kyx [N/m]', 'FontSize', 14);
% xlabel('N [Hz]','FontSize',14)
% ylabel('kyx [N/m]','FontSize',14)
% subplot(2, 2, 4);
% plot(N_vector_new, kyy_new, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
% title('kyy [N/m]', 'FontSize', 14);
% xlabel('N [Hz]','FontSize',14)
% ylabel('kyy [N/m]','FontSize',14)



% Plotting K
figure;
% Interpolation settings
N_vector_interp = linspace(min(N_vector), max(N_vector), 100); % Interpolation points
% Interpolated vectors
kxx_interp = interp1(N_vector, kxx, N_vector_interp, 'spline');
kxy_interp = interp1(N_vector, kxy, N_vector_interp, 'spline');
kyx_interp = interp1(N_vector, kyx, N_vector_interp, 'spline');
kyy_interp = interp1(N_vector, kyy, N_vector_interp, 'spline');
% Plot all graphs in one image
plot(N_vector, kxx, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
hold on;
plot(N_vector_interp, kxx_interp, 'r', 'LineWidth', 1.5);
plot(N_vector, kxy, '*-g', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'g');
plot(N_vector_interp, kxy_interp, 'm', 'LineWidth', 1.5);
plot(N_vector, kyx, '*-k', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'k');
plot(N_vector_interp, kyx_interp, 'c', 'LineWidth', 1.5);
plot(N_vector, kyy, '*-r', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'r');
plot(N_vector_interp, kyy_interp, 'b', 'LineWidth', 1.5);
hold off;
title('Stiffness [N/m] vs. speed [Hz]', 'FontSize', 14);
xlabel('N [Hz]','FontSize',14)
ylabel('kij [N/m]','FontSize',14)
legend('kxx', 'kxx interp', 'kxy', 'kxy interp', 'kyx', 'kyx interp', 'kyy', 'kyy interp', 'FontSize', 12, 'Location', 'best')



% Interpolation settings
N_vector_interp = linspace(min(N_vector), max(N_vector), 100); % Interpolation points

% Interpolated vectors
bxx_interp = interp1(N_vector, dxx, N_vector_interp, 'spline');
bxy_interp = interp1(N_vector, dxy, N_vector_interp, 'spline');
byx_interp = interp1(N_vector, dxy, N_vector_interp, 'spline');
byy_interp = interp1(N_vector, dyy, N_vector_interp, 'spline');

% Plotting Dxx
figure;
plot(N_vector, dxx, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
hold on;
plot(N_vector_interp, bxx_interp, 'r', 'LineWidth', 1.5);
hold off;

% Plotting Dxy
hold on;
plot(N_vector, dxy, '*-g', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'g');
plot(N_vector_interp, bxy_interp, 'm', 'LineWidth', 1.5);
hold off;

% Plotting Dyx
hold on;
plot(N_vector, dxy, '*-c', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'c');
plot(N_vector_interp, byx_interp, 'k', 'LineWidth', 1.5);
hold off;

% Plotting Dyy
hold on;
plot(N_vector, dyy, '*-r', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'r');
plot(N_vector_interp, byy_interp, 'b', 'LineWidth', 1.5);
hold off;

% Setting title and labels
title('Damping [N/(m/s)] vs. speed [Hz]', 'FontSize', 14);
xlabel('N [Hz]','FontSize',14)
ylabel('dij [N/(m/s)]','FontSize',14)

% Legend
legend('Dxx', 'Dxx interp', 'Dxy', 'Dxy interp', ...
       'Dyx', 'Dyx interp', 'Dyy', 'Dyy interp','FontSize', 12, 'Location', 'best')



% Plotting K
figure;
% Interpolation settings
N_vector_interp = linspace(min(N_vector), max(N_vector), 100); % Interpolation points
% Interpolated vectors
kxx_interp = interp1(N_vector, kxx, N_vector_interp, 'spline');
kxy_interp = interp1(N_vector, kxy, N_vector_interp, 'spline');
kyx_interp = interp1(N_vector, kyx, N_vector_interp, 'spline');
kyy_interp = interp1(N_vector, kyy, N_vector_interp, 'spline');
subplot(2, 2, 1);
plot(N_vector, kxx, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
hold on;
plot(N_vector_interp, kxx_interp, 'r', 'LineWidth', 1.5);
hold off;
title('kxx [N/m]', 'FontSize', 14);
xlabel('N [Hz]','FontSize',14)
ylabel('kxx [N/m]','FontSize',14)
subplot(2, 2, 2);
plot(N_vector, kxy, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
hold on;
plot(N_vector_interp, kxy_interp, 'r', 'LineWidth', 1.5);
hold off;
title('kxy [N/m]', 'FontSize', 14);
xlabel('N [Hz]','FontSize',14)
ylabel('kxy [N/m]','FontSize',14)
subplot(2, 2, 3);
plot(N_vector, kyx, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
hold on;
plot(N_vector_interp, kyx_interp, 'r', 'LineWidth', 1.5);
hold off;
title('kyx [N/m]', 'FontSize', 14);
xlabel('N [Hz]','FontSize',14)
ylabel('kyx [N/m]','FontSize',14)
subplot(2, 2, 4);
plot(N_vector, kyy, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
hold on;
plot(N_vector_interp, kyy_interp, 'r', 'LineWidth', 1.5);
hold off;
title('kyy [N/m]', 'FontSize', 14);
xlabel('N [Hz]','FontSize',14)
ylabel('kyy [N/m]','FontSize',14)



% Plotting D
figure;
% Interpolation settings
N_vector_interp = linspace(min(N_vector), max(N_vector), 100); % Interpolation points
% Interpolated vectors
bxx_interp = interp1(N_vector, dxx, N_vector_interp, 'spline');
bxy_interp = interp1(N_vector, dxy, N_vector_interp, 'spline');
byx_interp = interp1(N_vector, dxy, N_vector_interp, 'spline');
byy_interp = interp1(N_vector, dyy, N_vector_interp, 'spline');
subplot(2, 2, 1);
plot(N_vector, dxx, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
hold on;
plot(N_vector_interp, bxx_interp, 'r', 'LineWidth', 1.5);
hold off;
title('dxx [N/(m/s)]', 'FontSize', 14);
xlabel('N [Hz]','FontSize',14)
ylabel('dxx [N/(m/s)]','FontSize',14)
subplot(2, 2, 2);
plot(N_vector, dxy, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
hold on;
plot(N_vector_interp, bxy_interp, 'r', 'LineWidth', 1.5);
hold off;
title('dxy [N/(m/s)]', 'FontSize', 14);
xlabel('N [Hz]','FontSize',14)
ylabel('dxy [N/(m/s)]','FontSize',14)
subplot(2, 2, 3);
plot(N_vector, dxy, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
hold on;
plot(N_vector_interp, byx_interp, 'r', 'LineWidth', 1.5);
hold off;
title('dyx [N/(m/s)]', 'FontSize', 14);
xlabel('N [Hz]','FontSize',14)
ylabel('dyx [N/(m/s)]','FontSize',14)
subplot(2, 2, 4);
plot(N_vector, dyy, '*-b', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', 'b');
hold on;
plot(N_vector_interp, byy_interp, 'r', 'LineWidth', 1.5);
hold off;
title('dyy [N/(m/s)]', 'FontSize', 14);
xlabel('N [Hz]','FontSize',14)
ylabel('dyy [N/(m/s)]','FontSize',14)


