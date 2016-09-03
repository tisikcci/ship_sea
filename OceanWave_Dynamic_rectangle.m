clear all;
close all;
clc;
U =10;     %Wind_Speed
Afa = 8.1e-3;    %Constant
Beita = 0.74;     %Constant 
g = 9.8;             %Gravitational constant.
Omiga = 0.0001 : 0.1 : 2;
Omiga = Omiga.';
S = ((Afa*g^2)./(Omiga.^5)).*exp(-Beita*((g/U)./Omiga).^4); % P-M Spectrum
figure(1);
subplot(2,1,1),plot(Omiga./(2*pi),S);
subplot(2,1,2),plot(Omiga,S);
Omiga_Max = 8.565/U;
p = 0.5 + 0.82.*exp(-0.5.*((Omiga./Omiga_Max).^4));
q = 0.32.*exp(-0.5.*((Omiga./Omiga_Max).^4));
theata = -5e-4:2e-5:5e-4;
Extern = ones(1,length(theata));
G = (1/pi)*(1 + p*cos(2.*theata)+q*cos(4.*theata));  %SWOP—Stereo Wave Observation Project
S_extern = S*Extern; 
S_FW = S_extern.*G;
a = sqrt((2*0.5*0.5)*S_FW);    %Wave Amplitude

Omiga_extern = Omiga * Extern;
Initial_Phase = 2.*rand(length(Omiga),length(theata)); % 波浪起始点各种谐波的初相位

kn = Omiga.^2/g;  %波数

Wave_speed = (sqrt(g./kn)).';%不同谐波的波速（与波的频谱有关）

%仿真覆盖1km*1km区域的海面，如果要使海面都出现各次谐波的波涛起伏，
%海面区域以极坐标的形式描述，1000m半径，pi/2角度
%最波速最慢的谐波也必须传播到最远的区域;
t_start = 1000/ min(Wave_speed);

%海面点的描述
Sea = (10000e3:1:10001e3).'*theata;
R = (10000e3:1:10001e3).';

%计算某一时刻t的海面形状（t>t_Start）
t = t_start + 100;

%计算某一距离某一谐波传输所需要的时间
t_expand = (R*(1./Wave_speed)).';
%等效波源起始点振动时间
t_equal = t - t_expand;%(竖直方向表述距离，水平方向表述频率)

Z =  ones(size(Sea));
X = R*cos(theata);
Y = R*sin(theata);


for j = 1:length(theata);
    for i = 1:length(R);
        Z(i,j) = sum(a(:,j).*cos(Omiga.*t_equal(:,i) + Initial_Phase(:,j)));
    end;
end;
    
figure(2);
% pp = mesh(X,Y,Z);
pp = surfl(X,Y,Z,[-90 30],[0.5 0.6 1.5 10]);
colormap(cool);
shading interp;
hidden off;


axis([10000e3,10001e3,-5.0000e+003,5.0000e+003,-12,12]);

view([30,30,170]);

for time = 0:0.5:100;
    t_equal = t+time - t_expand;%(竖直方向表述距离，水平方向表述频率)
    for j = 1:length(theata);
        for i = 1:length(R);
            Z(i,j) = sum(a(:,j).*cos(Omiga.*t_equal(:,i) + Initial_Phase(:,j)));
        end;
    end;
    
    set(pp,'zdata',Z);
    drawnow;
end

    

        




