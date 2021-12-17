clear all
close all
%*******************************************************************************
%Script d'un signal en arche à partir d'un signal harmonique fondamental
%INPUT: A = amplitude
%INPUT: fe = fréquence d'échantillonnage
%INPUT: f = fréquence fondamentale
%*******************************************************************************

fe=100;
pas=1/fe;
A = 1;
f0 = 5;
%declaration axe du temps
t = 0:pas:1;
%decaleration axes des volts
y = A*sin(2*pi*f0*t);

%SIGNAL FONDAMENTAL
figure(1)
%declaration graphe
plot(t,y,'b')
title('Signal')
xlabel('temps en (s)')
ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)'])


%matrice harmoniques
H=48;
Mr = [];
for i=1:H
  ytemp = A*cos(2*pi*i*f0*t);
  Mr = [Mr;ytemp];
end

%matrice amplitudes
V=[];
for i=1:H
  Vtemp = -1/(i^2);
  V = [V;Vtemp];
end

%multiplication A et Mr
RR =[];
for i = 1:length(V)
  rtemp=V(i)*Mr(i,:);
  RR = [RR;rtemp];
end

%REPRESENTATION TRIDIMENSIONNELLE
figure(2)
M2 = cumsum(RR)
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
zlabel(['y = ',num2str(A),'*sin(2*pi*i*',num2str(f0),'*t)'])
ylabel('Nombre d''harmoniques')

h=1:H;

%DRECROISSANCE DES HARMONIQUES
figure(5)
stem(h,abs(V))
title('Décroissance harmoniques')
xlabel('Nombre d''harmoniques ')
ylabel('Amplitude')

%SOMME CUMULEE DES HARMONIQUES
figure(6)
plot(t, M2,'b')
title('Somme cumulée des harmoniques')
xlabel('temps en (s)')
ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)'])

%SIGNAL RECONSTITUE
figure(7)
plot(M2(H,:))
title('Signal reconstitué')
xlabel('temps en (s)')
ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)'])

%BILAN DE L ENSEMBLE DES REPRESENTATIONS
figure(8)
subplot(4,1,1), plot(t,y,'b') %signal harmonique fondamental
title('Signal')
AZ = title('Signal');
set(AZ,'fontsize',13)
xlabel('temps en (s)')
ZA=xlabel('temps en (s)')
set(ZA,'fontsize',7,'position',[0.5000 -1.3 -1])
ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)'])
AA=ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)']);
set(AA,'fontsize',7,'position',[-0.05 9.5367e-07 -1])
subplot(4,1,2), plot(t,RR,'b') %somme des harmoniques
title('Somme des harmoniques')
BZ = title('Somme des harmoniques');
set(BZ,'fontsize',13)
xlabel('temps en (s)')
ZB=xlabel('temps en (s)')
set(ZB,'fontsize',7,'position',[0.5000 -1.3 -1])
ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)'])
AB=ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)']);
set(AB,'fontsize',7,'position',[-0.05 6.0713e-07 -1])
subplot(4,1,3), plot(t, M2,'b') %somme cumulée des harmoniques
title('Somme cumulée des harmoniques')
CZ = title('Somme cumulée des harmoniques');
set(CZ,'fontsize',13)
xlabel('temps en (s)')
ZC=xlabel('temps en (s)')
set(ZC,'fontsize',7,'position',[0.5000 -2.45 -1])
ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)'])
AC=ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)']);
set(AC,'fontsize',7,'position',[-0.05 6.0713e-07 -1])
subplot(4,1,4), plot(M2(H,:)) %signal reconstitué
title('Signal reconstitué')
DZ = title('Signal reconstitué');
set(DZ,'fontsize',13)
xlabel('temps en (s)')
ZD=xlabel('temps en (s)')
set(ZD,'fontsize',7,'position',[60.0001 -2.45 -1.0000])
ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)'])
AD=ylabel(['y = ',num2str(A),'*sin(2*pi*',num2str(f0),'*t)']);
set(AD,'fontsize',7,'position',[-5.8 4.9251e-07 -1.0000])

%*******************************************************************************
%Script de la transformée de Fourier d'un signal harmonique
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
xlabel(['translation autour de la fréquence fondamentale (f0=',num2str(f0),')'])
ylabel(['Amplitude multipliée par la fréquence d''échantillonnage(fe=',num2str(fe),')'])

%LIEN ENTRE LA FFT ET LA SF
figure(10)
v=abs(V);
subplot(2,1,1), stem(h,v)
title('Décroissance harmoniques')
xlabel('Nombre d''harmoniques ')
ylabel('Amplitude')
subplot(2,1,2), stem(vec,2*abs(TF))
title('FFT')
xlabel(['translation autour de la fréquence fondamentale (f0=',num2str(f0),')'])
ylabel(['Amplitude * par la fréquence d''échantillonnage(fe=',num2str(fe),')'])
RT = ylabel(['Amplitude * par la fréquence d''échantillonnage(fe=',num2str(fe),')']);
set(RT,'fontsize',7)