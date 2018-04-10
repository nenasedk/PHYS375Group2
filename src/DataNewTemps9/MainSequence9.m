name = ['MSStar_', num2str(29),'.txt'];
fid = fopen(name);
M = cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f',[inf, 20], 'Delimiter', ',','HeaderLines',6));
%Mtotal = zeros(length(M(:,1)), length(M(1,:)),16 );
start = 1;
stop = 100;
X = zeros(stop,20);
Y = zeros(stop,1);
EmpiricLum = zeros(stop, 1);
empiricrad = zeros(stop, 1);
T_S = zeros(stop, 1);
Lsun = 3.828e26;
Rsun = 6.957e8;
sig = 5.67e-8;
Test1 = zeros(stop, 1);
Test2 = zeros(stop, 1);
Msun = 1.989e30;

for i = start:stop

name = ['MSStar_', num2str(i),'.txt'];
fid = fopen(name);
M = cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f',[inf, 20], 'Delimiter', ',','HeaderLines',6));
fclose(fid);

Mmax = max(M);
GAH = M(:,3);
Mmin = min(GAH (GAH>0));

X(i,:) = Mmax;

Y(i,:) = Mmin;
Test1(i,:) = Y(i,1);


    
T_S(i,:) = (X(i,5)./(4.*pi.*((X(i,1)).^2).*sig)).^(1/4);
 

end

Magnitude = 4.83 - 2.5.* log10((X(:,5)./Lsun));
SurfaceTemp = log10(Y);
SurfaceTemp1 = log10(T_S);
figure(1)
hold on
scatter(SurfaceTemp1, Magnitude)
set( gca, 'YDir', 'reverse' )
set(gca, 'XDir','reverse')
title('HR Diagram')
xlabel('Surface Temperature (logarithmic)')
ylabel('Absolute Magnitude')
saveas(gcf, 'HRdiagram.png')
hold off
%axis([160 20000 -10 10]) %reset the limits

Mass = X(:,4);
scMass = Mass./Msun;
scLum = X(:,5)./Lsun;
empLLsun1 = 0.35.*((Mass./Msun).^(2.62)); %M < 0.7Msun
empLLsun2 = 1.02.*((Mass/Msun).^(3.92)); % M > 0.7Msun
for i = start:stop
    if (Mass(i) < 0.7*Msun)
        EmpiricLum(i) = empLLsun1(i);
    
    elseif (Mass(i) > 0.7*Msun)
        EmpiricLum(i) = empLLsun2(i);
    end
end

figure(2)
hold on
logscLum = log10(scLum);
logscMass = log10(scMass);
scatter(logscMass, logscLum, 'DisplayName', 'Calculated Relation')
ylabel('L (Solar Units)')
xlabel('M (Solar Units')
title('Stellar Luminosity vs. Mass (Logarithmic Scale)')
logEmpiricLum = log10(EmpiricLum);
plot(logscMass, logEmpiricLum, 'DisplayName', 'Empirical Relation')
legend('Location','Northwest')
saveas(gcf, 'lummass.png')
hold off

%Calculating Radii
Radius = X(:,1)./Rsun;
logRad = log10(Radius);
for i = start:stop
    if (Mass(i) < 1.66*Msun)
        empiricrad(i) = 1.06.*((Mass(i)./Msun).^0.945);
        
    elseif (Mass(i) > 1.66*Msun)
        empiricrad(i) = 1.33.*((Mass(i)./Msun).^0.555);
    end
end

logempiricrad = log10(empiricrad);

figure(3)
hold on
scatter(logscMass, logRad, 'DisplayName', 'Calculated Relation')
plot(logscMass, logempiricrad, 'DisplayName', 'Empirical Relation')
xlabel('M (solar units)')
ylabel('R (solar units)')
legend('Location','Northwest')
title('Stellar Radius vs. Mass (Logarithmic Scale)')
saveas(gcf, 'radiusmass.png')
hold off

V = 1:100;

figure(4)
scatter(V, SurfaceTemp1)
title('Temp')

figure(5)
scatter(V, Magnitude)
title('magnitude')


