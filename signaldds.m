clear all
close all
%*******************************************************************************
%Script d'un signal en dents de scie � partir d'un signal harmonique fondamental
%INPUT: A = amplitude
%INPUT: fe = fr�quence d'�chantillonnage
%INPUT: f = fr�quence fondamentale
%*******************************************************************************
%SIGNAL FONDAMENTAL
fe=100;     %Fr�quence �chantillonnage
A = 1;      %Amplitude du signal
f0 = 5;     %Fr�quence
pas=1/fe;
%Declaration axe du temps
t = 0:pas:1;
%Declaration axes des volts
y = A*sin(2*pi*f0*t);

figure(1)
%declaration graphe
plot(t,y,'b')
%declaration titre et axes
title('Signal')
xlabel('temps en (s)')
ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)'])

%matrice harmoniques
H=48;
Mr = [];
for i=1:H
  ytemp = A*sin(2*pi*i*f0*t);
  Mr = [Mr;ytemp];
end

%Vecteurs des coeeficients k
V=[];
for i=1:H
  Vtemp = (-2*A/(i*pi));
  V = [V;Vtemp];
end

%multiplication  des coefficients A et de la matrice des harmoniques
RR =[];
for i = 1:length(V)
  rtemp=V(i)*Mr(i,:);
  RR = [RR;rtemp];
end

%REPRESENTATION TRIDIMENSIONNELLE 
figure(2)
%Somme cumul�e
M2 = cumsum(RR);
%tracer graphe 3d
mesh(M2)
title(['Signal en fonction du nombre d''harmoniques et du temps'])
xlabel('temps en (s)')
zlabel(['z = ',num2str(A),'*sin(2*pi*i*',num2str(f0),'*t)'])
ylabel('Nombre d''harmoniques')

%SOMME DES HARMONIQUES
figure(3)
plot(t,RR,'b')
title('Somme des harmoniques')
xlabel('temps en (s)')
ylabel(['y = ',num2str(A),'*sin(2*pi*i*',num2str(f0),'*t)'])

%REPRESENTATION TRIDIMENSIONNELLE AMPLITUDE TEMPS NOMBRE D HARMONIQUES
figure(4)
mesh(RR)
title(['Signal en fonction du nombre d''harmoniques et du temps'])
xlabel('temps en (s)')
zlabel(['z = ',num2str(A),'*sin(2*pi*i*',num2str(f0),'*t)'])
ylabel('Nombre d''harmoniques')

%DRECROISSANCE DES HARMONIQUES
h= 1:H;
figure(5)
stem(h,abs(V))
title('D�croissance harmoniques')
xlabel('Nombre d''harmoniques ')
ylabel('Amplitude')

%SOMME CUMULEE DES HARMONIQUES
figure(6)
plot(t, M2,'b')
title('Somme cumul�e des harmoniques')
xlabel('temps en (s)')
ylabel(['y = ',num2str(A),'*sin(2*pi*i*',num2str(f0),'*t)'])

%SIGNAL RECONSTITUE
figure(7)
plot(M2(H,:))
title('Signal reconstitu�')
xlabel('temps en (s)')
ylabel(['y = ',num2str(A),'*sin(2*pi*i*',num2str(f0),'*t)'])

%BILAN DE L ENSEMBLE DES REPRESENTATIONS
figure(8)
subplot(4,1,1), plot(t,y,'b') %signal harmonique fondamental
title('Signal')
AZ = title('Signal');
set(AZ,'fontsize',13)
xlabel('temps en (s)')
ZA=xlabel('temps en (s)')
set(ZA,'fontsize',7,'position',[0.5000 -5 -1])
ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)'])
AA=ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)']);
set(AA,'fontsize',7,'position',[-0.05 9.5367e-07 -1])
subplot(4,1,2), plot(t,RR,'b') %somme des harmoniques
title('Somme des harmoniques')
BZ = title('Somme des harmoniques');
set(BZ,'fontsize',13)
xlabel('temps en (s)')
ZB=xlabel('temps en (s)')
set(ZB,'fontsize',7,'position',[0.5000 -5 -1])
ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)'])
AB=ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)']);
set(AB,'fontsize',7,'position',[-0.05 6.0713e-07 -1])
subplot(4,1,3), plot(t, M2,'b') %somme cumul�e des harmoniques
title('Somme cumul�e des harmoniques')
CZ = title('Somme cumul�e des harmoniques');
set(CZ,'fontsize',13)
xlabel('temps en (s)')
ZC=xlabel('temps en (s)')
set(ZC,'fontsize',7,'position',[0.5000 -7 -1])
ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)'])
AC=ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)']);
set(AC,'fontsize',7,'position',[-0.05 6.0713e-07 -1])
subplot(4,1,4), plot(M2(H,:)) %signal reconstitu�
title('Signal reconstitu�')
DZ = title('Signal reconstitu�');
set(DZ,'fontsize',13)
xlabel('temps en (s)')
ZD=xlabel('temps en (s)')
set(ZD,'fontsize',7,'position',[60.0001 -6.5 -1.0000])
ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)'])
AD=ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)']);
set(AD,'fontsize',7,'position',[-5.8 4.9251e-07 -1.0000])

%*******************************************************************************
%Script de la transform�e de Fourier d'un signal harmonique
% FFT algo qui permet de retrouver la TF d'un signal
% Signal : TF = fft(signal);
%*******************************************************************************


%FFT
signal = M2(H,:);
tf = fft(signal);
Vect_freq = 0:1:fe;
f=floor(fe/2);
vec=Vect_freq(1,1:f);
TF=tf(1,1:floor(fe/2));
figure(9)
stem(vec,2*abs(TF))
title('FFT')
xlabel(['translation autour de la fr�quence fondamentale (f0=',num2str(f0),')'])
ylabel(['Amplitude multipli�e par la fr�quence d''�chantillonnage(fe=',num2str(fe),')'])

figure(10)
v=abs(V);
subplot(2,1,1), stem(h,v)
title('D�croissance harmoniques')
xlabel('Nombre d''harmoniques ')
ylabel('Amplitude')
subplot(2,1,2), stem(vec,2*abs(TF))
title('FFT')
xlabel(['translation autour de la fr�quence fondamentale (f0=',num2str(f0),')'])
ylabel(['Amplitude * par la fr�quence d''�chantillonnage(fe=',num2str(fe),')'])
RT = ylabel(['Amplitude * par la fr�quence d''�chantillonnage(fe=',num2str(fe),')'])
set(RT,'fontsize',7)