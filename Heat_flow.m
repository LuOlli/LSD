clear all
load('TvsD.mat')%(TvsD(:,1) = temperature, TvsD(:,2) = depth )(degs C, meters)
load('CvsD.mat')%(CvsD(:,1) = conductivity, CvsD(:,2) = depth )(W/(m⋅K), meters)
%a)
%Calculate heat flow for layers of 100 m up to 900 m
layer_width = 100;%m
Depth = zeros(9,1);
DeltaT = zeros(9,1);
ThCond = zeros(9,1);
Heatflow = zeros(9,1);
gradT = zeros(9,1);
for n=1:9
    clear temp_TvsD temp_CvsD
    lb = n*layer_width-layer_width;
    ub = n*layer_width;
    Depth(n,:)=ub;
    temp_TvsD=TvsD((TvsD(:,2)> lb)&(TvsD(:,2)<= ub),:);%values between lb < D <= ub
    temp_CvsD=CvsD((CvsD(:,2)> lb)&(CvsD(:,2)<= ub),:);
    DeltaT(n,:)=temp_TvsD(end,1)-temp_TvsD(1,1);%temperature difference in the layer
    ThCond(n,:) = mean(temp_CvsD(:,1));%Average of the measured conductivities for each step
    gradT(n,:)  = DeltaT(n,:)/layer_width;
    Heatflow(n,:) = gradT(n,:)*ThCond(n,:);% heat flow for a step
end
disp('   Delta T|Conductivity|Heatflow')
disp([DeltaT ThCond Heatflow])
figure(1)
hold on
plot(abs(Heatflow),Depth,'ko')
set(gca, 'YDir','reverse')
ylabel('Depth(m)')
xlabel('Heat flow (W m^{-2})')
ylim([0,1000])

%b)

% Depth = linspace(0,15000)';
DeltaT2 = 5;%C
s=1.5e-6;%m2s-1
t2=11e3*365*24*60*60;
t1=100e3*365*24*60*60;
DeltaTzt = DeltaT2.*(erfc(Depth./sqrt(4*s*t2))-erfc(Depth./sqrt(4*s*t1)));
Deltagzt = -DeltaT2.*((pi*s*t1)^(-0.5)*exp(-Depth.^2/(4*s*t1))-(pi*s*t2)^(-0.5)*exp(-Depth.^2/(4*s*t2)));
% plot(Depth,DeltaTzt,'k')
% plot(Depth,Deltagzt,'k')
% xlabel('Depth(m)')
% ylabel('Gradient change (mK/m)')
Heatflow_cor  = ThCond.*Deltagzt+Heatflow
plot(Heatflow_cor,Depth,'ro')
hold off
legend('Heat flow','Corrected heat flow')
ylim([0,1000])
%% 2
z = linspace(0,1e4)';
Q0 = Heatflow_cor(1);
Lambda = mean(CvsD(:,1));
T0=4;%C
A=3e-6;%Wm^-3
T = T0 + (Q0/Lambda).*z-(A/(2*Lambda)).*z.^2;
figure(2)
hold on
plot(T,z,'k-')
set(gca, 'YDir','reverse')
ylabel('Depth(m)')
xlabel('Temperature(C^\circ)')
%b)
Q0b = [Q0-10e-3 Q0 Q0+10e-3];%Varied heat flow +-10mWm^-2
Lambdab = [Lambda-1 Lambda Lambda+1];%Varied thermal conductivity +-1Wm^1K^-1
Ab =  [A-1e-6 A A+1e-6];
nn = 1;
figure(3)
hold on
for i=1:3
    for j=1:3
        for k=1:3
            Tb(:,nn) = T0 + (Q0b(i)/Lambdab(j)).*z-(Ab(k)/(2*Lambdab(j))).*z.^2;
            plot(Tb(:,nn),z,'k')
            nn=nn+1;
        end
    end
end
set(gca, 'YDir','reverse')
ylabel('Depth(m)')
xlabel('Temperature(C^\circ)')

%% 3
clear all
L = 50e3;%m
TL = 500;%C
u = [-1e4 -1e3 -1e2 -1e1 -1 1 1e1 1e2 1e3 1e4]'/((1e+6)*360*24*60*60)';%m/s
z = linspace(0,10e3);
s=1.5e-6;%m2s-1, I will use this value from 1.
T=(TL*(1-exp(-u*z/s))./(1-exp(-u*L/s)))';
plot(T,z)
set(gca, 'YDir','reverse')
% set(gca, 'XScale', 'log')
ylabel('Depth(m)')
xlabel('Temperature(C^\circ)')
%b)
ThCond = 2.5;
Tb = T(:,6);
grad = diff(Tb)./diff(z');






