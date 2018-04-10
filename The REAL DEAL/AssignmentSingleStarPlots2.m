Msun = 1.99e30;
G = 6.67e-11;
Lsun = 3.828e28;

name = 'MSStar_88.txt';
fid = fopen(name);
M = cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f',[inf, 20], 'Delimiter', ',', 'HeaderLines',11));
%Mtotal = zeros(length(M(:,1)), length(M(1,:)),3 );
fid = fopen(name);
fclose(fid);

Rmax = max(M(:,1)); %Highest radius given by the code
Dc = max(M(:, 2)); %Highest density given by the code
Tx = M(:,3); %All the temperatures all the way through the star
Tc = max(M(:,3)); %Core temp found by picking out the highest T 
Mmax = max(M(:, 4)); %Total Mass
Lmax = max(M(:, 5)); %Max Luminosity
R = M(:,1)./Rmax; %Normalized Radius
Dens = M(:, 2)./Dc; %Normalized Density
T = Tx./Tc; %Normalized Temperature
Mass = M(:,4)./Mmax; %Normalized Mass
L = M(:,5)./Lmax; %Normalized Luminosity
OD = log10(M(:,6)); %Optical depth if you're interested
dLdr = M(:,17); %dL/dr is given by the code

TotP = M(:,16); %Total P, also given by the code
Pmax = max(TotP); %Maximum Pressure
P = TotP./Pmax; %Normalized Pressure

dTrad = M(:,19); %These two are exactly what they say, output by the code
dTcon = M(:,20);

dTdr = min(dTrad, dTcon); %Pick the min between the two to represent dTdr

dPdr = -((G.*(M(:,4)).*(M(:,2)))./(M(:,1)));
dPdT = M(:,18);

PTplot = ((Tx)./(TotP)).*dPdr.*((dTdr).^(-1)); %Math version of dlogP/dlogT


Deg = M(:,13)./Pmax; %Normalized by maximum pressure degeneracy, gas and radiation
Gas = M(:,14)./Pmax;
Rad = M(:,15)./Pmax;

k1 = log10(M(:,7)); %kff
k2 = log10(M(:,8)); %kH
k3 = log10(M(:,9));%kes

%effective opacity
keff = log10(((M(:,8)).^(-1) + (max(M(:,7),M(:,9))).^(-1)).^-1); %the effective opacity


%Calculate dlogP and dlogT

logP = log10(M(:,16));
logT = log10(Tx);
dlogP = zeros(length(logP), 1);
dlogT = zeros(length(logT),1);

for i = 1:(length(logP)-1) 
   dlogP(i) = (logP(i+1)) - (logP(i));
   dlogT(i) = logT(i+1) - logT(i); 
end
%then to plot you want to plot the following
dlogPdlogT = dlogP./dlogT;

CalculatedSurfaceTemperature = (Lmax/(4*pi*(5.67e-8)*(Rmax^2)))^(1/4);
%PLOTS

figure(1) %the four-line plot, normalized
plot(R,Dens, 'DisplayName', '\rho/\rho_c')
hold on
title('Internal Normalized Star Structure of a High Mass Star (M = 2.5M_s_u_n)')
plot(R, T,'DisplayName', 'T/T_c')
plot(R, Mass, 'DisplayName', 'M/M_T_o_t')
plot(R,L, 'DisplayName', 'L/L_m_a_x' );
%plot(R, OD, 'DisplayName', 'Opticial Depth')
%plot(R, P, 'DisplayName', 'Pressure')
xlabel('Radius (R/Rmax)')
ylabel('\rho/\rho_c, T/T_c, M/M., L/L.')
legend()
saveas(gcf, 'starstructure.png')
hold off


figure(2) %Temperature
hold on
plot(R,M(:,3))
ylabel('Temperature (K)')
xlabel('Radius (R/Rmax)')
title('Temperature vs. Radius')
saveas(gcf, 'temp.png')
hold off

figure(3) %Mass
hold on
plot(R,M(:,4))
ylabel('Mass (kg)')
xlabel('Radius (R/Rmax)')
title('Mass vs. Radius')
saveas(gcf, 'mass.png')
hold off

figure(4) %Luminosity
hold on
plot(R,M(:,5))
title('Luminosity vs. Radius')
ylabel('Luminosity (W)')
xlabel('Radius (R/Rmax)')
saveas(gcf, 'Luminosity.png')
hold off

figure(5) %Density
hold on
plot(R,M(:,2))
title('Density vs. Radius')
ylabel('Density (kg/m^{3})')
xlabel('Radius (R/Rmax)')
saveas(gcf, 'dens.png')
hold off

figure(6) %Plots Normalized pressures
hold on
plot(R, P, 'DisplayName', 'Total Pressure')
plot(R, Deg,'DisplayName', 'Degeneracy Pressure')
plot(R, Gas,'DisplayName', 'Gas Pressure')
plot(R, Rad,'DisplayName', 'Radiation Pressure')
xlabel('R/Rmax')
ylabel('Pressure')
title('Normalized Pressures vs. Normalized Radius')
legend()
saveas(gcf, 'pressures4.png')
hold off

figure(7) %Plots total pressure (with actual values
hold on
plot(R, TotP, 'DisplayName', 'Total Pressure')
xlabel('R/Rmax')
ylabel('Pressure (Pa)')
title('Total Pressure vs. Radius')
saveas(gcf, 'totalpressure.png')
hold off

figure(8) %Plots dL/dr
hold on
plot(R,dLdr, 'DisplayName','dL/dr')
title('dL/dr vs. R/R_m_a_x')
xlabel('R/R_m_a_x')
ylabel('dL/dr (W/m)')
saveas(gcf, 'dLdr.png')
hold off

figure(9) %Plots Energy Generation Rates
hold on
plot(R, M(:,10), 'DisplayName','PP')
plot(R, M(:,11), 'DisplayName','CNO')
plot(R,M(:,12), 'DisplayName','3\alpha')
legend()
ylabel('\epsilon_P_P, \epsilon_C_N_O, \epsilon_3_\alpha (W/kg)')
xlabel('R/Rmax')
title('Energy Generation Rates vs. Radius')
saveas(gcf, 'EGR.png')
hold off

figure(10) %all the commented plots are the previous less smooth iterations
hold on
smoothed = smooth(dlogPdlogT);
doublesmoothed = smooth(smoothed);
supersmoothed = smooth(doublesmoothed);
%plot(R,smoothed, 'DisplayName', 'Smoothed')
%plot(R(1:50:end), dlogPdlogT(1:50:end), 'DisplayName', 'Unsmoothed')
%plot(R, doublesmoothed, 'DisplayName', 'Doublesmoothed')
plot(R,supersmoothed, 'DisplayName', 'Supersmoothed')
title('dlogP/dlogT vs. Radius')
ylabel('dlogP/dlogT (dlog(Pa)/dlog(K))')
xlabel('R/Rmax')
saveas(gcf, 'dlogPdlogT.png')
hold off

figure(11) %Opacities
hold on
plot(R, k1,'DisplayName', '\kappa_f_f')
plot(R, k2,'DisplayName', '\kappa_H')
plot(R, k3,'DisplayName', '\kappa_e_s')
plot(R, keff, 'DisplayName', '\kappa_e_f_f')
legend('location','northwest')
title('Log plot of Opacities vs. Radius')
ylabel('Opacities (log(m^{2}/kg))')
xlabel('R/Rmax')
ylim([-3 10])
saveas(gcf, 'opacity.png')
hold off
