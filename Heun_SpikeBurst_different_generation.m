% Program created by Bellante Luca and Giampieri Andrea starting from 
% the article by Professor Deng of the University of Nebraska.
% Univpm, Febraury 2022


close all 
clear all
clc
 
n = 30000;
h =  0.001;
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

 
%Dati sul circuito, presi dalla fig.8
Gna = 0.17;
Dna = -0.06;
i1 = 0.5;
i2 = 1;
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
e = 0.01;

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
Vna(I) = I/Gna + p/Dna * (((atan((I-im)/p))) + ((atan(im/p))));
%corrente canale potassio N-shape
Ik(V) = Gk*V+Dk*u*((atan((V-vm)/u))) + (Dk*u*((atan(vm/u))));


% Applico il sistema risolvente del eq.13 

% Pagina 22 in poi 
Ana(1)=Ga*Ena;
Inap(1)=-Ga*Ena+Iext;
%Vc(1) = (1/Gna)*Inap(1) + Ena;
Vc(1)=Ena*(Gna/(Gna+Ga))+Iext/(Gna+Ga);
Ak(1)=-Ga*Ek;
Ikp(1)=-Ga*Ek+Iext;
% Ak(1)=0; %la cond iniziale aumenta la frequenza delle spike
% Ana(1)=8;

t = 1;
for i=1 : n
    
     if(i==300*t)
        clc;
        disp("Simulation progress: " + string((i/n)*100) + "%");
        t = t+1;
     end 
    
    tempo(i+1)=tempo(i)+h;
    
% Condizione necessaria e sufficiente per l'interruzione dello Spike Burst
    if((i >=6000 && i <= 8000) || (i >= 15000 && i <= 17000) || (i >= 24000 && i <= 26000 ))
        % Questa condizione l'ho dedotta dal sistema risolvente, sulla base del fatto
        % che essa ha lo scopo di "Eccitare" il bipolo, ossia il neurone
        % artificiale, quindi va da se che il valore di essa dipende solo
        % ed esclusivamente da un ipotetico Gc controllato appunto da me
        Iext = Inap(i) + Ikp(i) + Ana(i) - Ak(i);
        % Impongo questa condizione poichè va da se che se vale la
        % precedente allora la tensione calcolata al istante Vc(i)
        % =-(Inap(i) + Ikp(i) + Ana(i) - Ak(i) - Iext)/C == 0 , Quindi 
        % la esplicito nella relazione.
        % Più precisamente, considerando la derivata si ha:
        % C*Vc(i+1)=0 <--> Vc(i)=0 (osservare che la relazione Iext vale nel 
        % istante i) 
         Vc(i) = 0;
    end 
    
    
%     % Corrente totale canale passivo
     Ip(i) = Ikp(i) + Inap(i);
%     
%     % Somma e differenza delle correnti della pompa ionica 
     Is(i) = Ana(i) + Ak(i);
     Ia(i) = Ana(i) - Ak(i);
    
    
     
    %---------------  INIZIO SISTEMA RISOLVENTE  -----------------------%
    
    
    % Applico il metodo di Heun alla Vc
    % Iext l'ho considerata come una costante nulla(=0)
    F1vc = -(Inap(i) + Ikp(i) + Ana(i)-Ak(i) - Iext)/C ;
    F2vc = -((Inap(i)+h*F1vc)+(Ikp(i)+h*F1vc)+(Ana(i)+h*F1vc)-(Ak(i)+h*F1vc) - Iext)/C; % NoN so se la Iext è corretta
    Vc(i+1)= Vc(i) + (h/2)*F1vc + (h/2)*F2vc;
    
    % Applico il metodo di Heun alla Ana
    F1ana = l*Ana(i)*(Vc(i)-y*(Ana(i)-Ak(i)));
    F2ana = l*(Ana(i)+h*F1ana)*((Vc(i)+h*F1ana)-y*((Ana(i)+h*F1ana)-(Ak(i)+h*F1ana)));
    Ana(i+1) = Ana(i) + (h/2)*F1ana + (h/2)*F2ana;
    
    % Applico il metodo di Heun alla Ak
    F1ak = l*Ak(i)*(-Vc(i) + y*(Ana(i)-Ak(i)));
    %controllare segno Vc , possibile errore con il meno
    F2ak = l*(Ak(i)+h*F1ak)*(-(Vc(i)+h*F1ak)+y*((Ana(i)+h*F1ak)-(Ak(i)+h*F1ak)));
    Ak(i+1) = Ak(i) + (h/2)*F1ak + (h/2)*F2ak;
    
   % Relazione costante, vedere pagina 16 del pdf di BoDeng
    Inap(i+1) = Inap(i);
  %  Inap(i+1)=Gna*(Vc(i+1)-Ena);
     
 %-------------------  FINE SISTEMA RISOLVENTE  -----------------------%
   
    
   
    
    Vnap(i+1)=Vna(Inap(i+1));
    Ikp(i+1)=Ik(Vc(i+1)-Ek);
end 


% Ip(n+1) = Ikp(n+1); % somma correnti canali passivi
% Ia(n+1)=Ana(n+1)-Ak(n+1); %corrente attiva pompa ionica
% Is(n+1)=Ana(n+1)+Ak(n+1); %somma correnti pompa ionica


plot(tempo,Vc);
title({'','Vc / time',''});

