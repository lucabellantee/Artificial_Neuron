% Program created by Bellante Luca and Giampieri Andrea starting from 
% the article by Professor Deng of the University of Nebraska.
% Univpm, Febraury 2022


clear all
close all
clc
 
n = 1400;
T = 0.01;
tempo = zeros(1,n);
Vc = zeros(1,n);
Ana = zeros(1,n);
Ak = zeros(1,n);
Ip = zeros(1,n);
Ikp = zeros(1,n);
Vnap = zeros(1,n);
Inap = zeros(1,n);
Ia = zeros(1,n);
Is = zeros(1,n);
Ic = zeros(1,n);
Ika = zeros(1,n);
Inaa = zeros(1,n);
 
%Dati sul circuito, presi dalla fig.10
Gna = 0.17;
Dna = -0.06;
i1 = 0.1;
i2 = 0.3;
Ena = 0.6;
Gk = 1;
Dk = -1.25;
v1 = 0.5;
v2 = 2;
Ek = -0.7;
C = 0.01;
l = 0.05;
y = 0.1; Ga=1/y;
Iext = 0;
e = 0.001;

% % Parametri con i quali confrontare la Iext
% Imthr=(Gna+Gk+Ga)*(Ek+v1)-(Gna*Ena+Gk*Ek)
% Inathr=(i1*(1+y*Gna)/(y*Gna))+(Ena/y)
% Ikthr=(v1*y/(1+Gk))+(Ena/y)
 
im=(i1+i2)/2;%corrente media canale sodio
p=sqrt(abs(Dna)/(Gna+Dna))*(i2-i1)/2; %parametro rho canale sodio
vm=(v1+v2)/2; %tensione media canale potassio
u=sqrt(Gk/abs(Gk+Dk))*(v2-v1)/2; %parametro mu canale potassio
syms Vna(I) Ik(V)

%tensione canale sodio S-shape
%Vna(I) = I/Gna + p/Dna * (((atan((I-im)/p))) + ((atan(im/p)))); %smooth
Vna(I) = piecewise((I<=i1),I/Gna,(i1<I) & (I<i2),I/Gna +(I-i1)/Dna,(I>=i2),I/Gna + (i2-i1)/Dna);

%corrente canale potassio N-shape
%Ik(V) = Gk*V+Dk*u*((atan((V-vm)/u))) + (Dk*u*((atan(vm/u)))); %smooth
Ik(V) = piecewise((V<=v1),Gk*V,(v1<V) & (V<v2),Gk*V +(V-v1)*Dk,V>=v2,Gk*V +(v2-v1)*Dk );

%corrente canale potassio N% Pagina 22 in poi, varia la freq dei picchi di vc
Vc(1)=((Gna*Ena) +(Gk*Ek))/(Gna+Gk+Ga) + Iext(1)/(Gna+Gk+Ga);
Ana(1)=2;
%Inap(1)=-Ga*Ena+Iext(1);
% Vc(1) = (1/Gna)*Inap(1) + Ena;
Ak(1)=Ana(1)+Vc(1)*Ga;
%Ikp(1)=-Ga*Ek+Iext(1);
%Ak(1)=8; %la cond iniziale aumenta la frequenza delle spike
%Ana(1)=8;
 

% Condizioni iniziali
Inaa(1)=Inap(1)+Ana(1); %somma dei canali sodio
Ika(1)=Ikp(1)+Ak(1); %somma dei canali potassio
Ia(1)=Ana(1)-Ak(1); %corrente attiva pompa ionica
Is(1)=Ana(1)+Ak(1); %somma correnti pompa ionica
Inap(1)=Gna*(Vc(1)-Ena); %corrente canale sodio pag 28
Vnap(1)=Vna(Inap(1));
Ikp(1)=Ik(Vc(1)-Ek);

t = 1;
for i=1:n
    
    if(i== 14*t)
        clc;
        disp("Simulation progress: " + string((i/n)*100) + "%");
        t = t+1;
    end
    
tempo(i+1)=tempo(i)+T;
 
%                       Sistema Risolvente:
Vc(i+1) = (Ak(i)*T - Ana(i)*T + C*Vc(i) + Iext*T - Inap(i)*T -Ikp(i)*T)/C;
Ana(i+1) = Ana(i)*(T*Vc(i)*l + Ak(i)*T*l*y - Ana(i)*T*l*y + 1);
Ak(i+1) = -Ak(i)*(T*Vc(i)*l + Ak(i)*T*l*y - Ana(i)*T*l*y - 1);

Inap(i+1)=Gna*(Vc(i+1)-Ena); %corrente canale sodio pag 28
Vnap(i+1)=Vna(Inap(i+1));
Ikp(i+1)=Ik(Vc(i+1)-Ek);

Ia(i+1)=Ana(i+1)-Ak(i+1); %corrente attiva pompa ionica
Is(i+1)=Ana(i+1)+Ak(i+1); %somma correnti pompa ionica
Inaa(i+1)=Inap(i+1)+Ana(i+1); %somma dei canali sodio
Ika(i+1)=Ikp(i+1)+Ak(i+1); %somma dei canali potassio

Vnap(i+1)=Vna(Inap(i+1));
Ikp(i+1)=Ik(Vc(i+1)-Ek);

Ia(i+1)=Ana(i+1)-Ak(i+1); %corrente attiva pompa ionica
Is(i+1)=Ana(i+1)+Ak(i+1); %somma correnti pompa ionica
Inaa(i+1)=Inap(i+1)+Ana(i+1); %somma dei canali sodio
Ika(i+1)=Ikp(i+1)+Ak(i+1); %somma dei canali potassio
end

 
subplot(3,4,1);
fplot (Ik); hold on;
title({'','Ik / Vk',''});
subplot(3,4,2);
fplot (Vna)
title({'','Vna / Ina',''});
subplot(3,4,3);
plot (tempo,Ikp)
title({'','Ikp / tempo',''});
subplot(3,4,4);
plot (tempo,Inap)
title({'','Inap / tempo',''});
subplot(3,4,5);
plot (tempo,Vnap)
title({'','Vnap / tempo',''});
subplot(3,4,6);
plot (tempo,Ana)
title({'','Ana / tempo',''});
subplot(3,4,7);
plot (tempo,Ak)
title({'','Ak / tempo',''});
subplot(3,4,8);
plot (tempo,Vc)
title({'','Vc / tempo',''});
subplot(3,4,9);
plot (tempo,Ia)
title({'','Ia / tempo',''});
subplot(3,4,10);
plot (tempo,Is)
title({'','Is / tempo',''});
subplot(3,4,11);
plot (tempo,Ika)
title({'','Ik all / tempo',''});
subplot(3,4,12);
plot (tempo,Inaa)
title({'','INa all / tempo',''});