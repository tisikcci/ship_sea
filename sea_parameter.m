%----------------------------海水参数计算---------------------------------%
%根据海水的温度和经验公式，计算海水的盐浓度、复介电常数、张弛时间、
%海水媒质的电导率、海水静态介电常数、

%
%%
T = 20; %海水温度
Ss = 32.54e-3; %海水盐浓度
if(Ss<0||Ss>40e-3)
    disp('请输入复合要求的海水浓烟度！');
end
%已知参数
emsn0 = 1/(36*pi)*1e-9;      %自由空间介电常数
a = 0;                       %张弛时间分布的经验数
emsnoo = 4.9;                %海水在无限大频率时的介电常数
emsn_T_0 = 87.134-1.949e-1*T-1.276e-2*T^2+2.491e-4*T^3; %纯水的介电常数随温度的变化
a_T_Ss = 1.0+1.613e-5*T*Ss-3.656e-3*Ss+3.21e-5*Ss^2-4.232e-7*Ss^3;
%tao_T_0 = 1/(2*pi)*(1.1109e-10-3.824e-1*T-1.276e-2*T^2+2.491e-4*T^3);%纯水的张弛时间
tao_T_0 = 1.768e-11-6.086e-13*T+1.104e-14*T^2-8.111e-17*T^3;
b_T_Ss = 1.0+2.282e-5*T*Ss-7.638e-4*Ss-7.76e-6*Ss^2+1.105e-8*Ss^3;
emsn_s_T_Ss = emsn_T_0*a_T_Ss;    %海水静态介电常数
tao_T_Ss = tao_T_0*b_T_Ss;        %张弛时间

deta = 25-T;
deta_25_Ss = Ss*(0.18252-1.4619e-3*Ss+2.093e-5*Ss^2-1.282e-7*Ss^3);%25摄氏度时海水的电导率
afa = 2.033e-2+1.266e-4*deta+2.464e-6*deta^2-...
          Ss*(1.849e-5-2.551e-7*deta+2.551e-8*deta^2);
deta_T_Ss = deta_25_Ss*exp(-deta*afa);         %海水电导率
%f=0.5e9:0.1e9:20e9;                     %频率采样
%f = 1e9:1e9:1000e9;
f = 5.5e9;
emsn = emsnoo+(emsn_s_T_Ss-emsnoo)./(1+(2j*pi*f*tao_T_Ss).^(1-a))-1j*deta_T_Ss./(2*pi*f*emsn0); %海水介电常数
real_emsn = real(emsn);
imag_emsn = abs(imag(emsn));
%imag_emsn1 = 2*pi*f*tao_T_Ss*(emsn_s_T_Ss-emsnoo)./(1+(2*pi*f*tao_T_Ss).^(1-a));
 %%
% figure(1)  
% plot(f./1e9,real_emsn)
% hold on
 
 %figure(2)
 %plot(f./1e9,imag_emsn)
 %hold on
 
 %figure(3)
 %plot(f./1e9,imag_emsn1)
 
   

