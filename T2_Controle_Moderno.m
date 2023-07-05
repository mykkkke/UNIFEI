pkg load control;
pkg load signal;
graphics_toolkit 'gnuplot'

clear all;
close all;
clc
%Parâmetros
m = 0.1;
k = 10;
g = 9.81;
b = 0.1;
h = 0.01;
F = 0;
ta = 0.1;

A = [ 0 1; -k/m -b/m];
B = [0; 1/m];
C = [1 0];
D = [0];

sys = ss(A, B, C, D);

% Obtenção da função de transferência
[num, den] = ss2tf(A, B, C, D);

% Criação do objeto de função de transferência
tf_sys = tf(num, den);
#bode(tf_sys);
#margin(tf_sys);

%Polos
p1 = -2.2857;
p2 = -40;
p3 = -60;

%Matrizes Expandidas
A_a = [A [0;0]; -C 0];
B_b = [B; 0];

kk = acker(A_a, B_b, [p1 p2 p3]); %Ganho Sistema
k1 = kk(1, 1:2); %Ganhos Retroação
jj = -kk(1, 3); %Ganhos Expensão


%Matrizes com Ganhos Aplicados
A_aa = [A-B*k1 B*jj; -C 0];
B_bb = [0;0;20];
C_cc = [C 0];

printf("  Matriz A\n");
disp(A_aa)
printf("\n");


t1 = ss(A_aa, B_bb, C_cc); % Malha Fechada do Sistema Expandido
[n, d] = ss2tf(A_aa, B_bb, C_cc, 0) % Passando para Função de Transferencia
tf_sy1 = tf(n, d) %Junta denom, numerador



L = acker(A', C', 1/3*[p1 p2]); %Ganho Observador Continuo

%Matrizes com Retroação de Estados Observados com uma Expansão de Polos
A_aaa = [(A) -B*k1 B*jj;L'*C A-B*k1-L'*C B*jj; -C 0 [0 0]];
B_bbb = [0; 0; 0; 0; 20];
C_ccc = [C 0 0 0];
t2 = ss(A_aaa, B_bbb, C_ccc,0)
[n, d] = ss2tf(t1);
tf_sy1 = tf(n, d)

figure;
step(t2);

t = [0:0.01:3];  % Vetor de tempos para calcular a resposta do sistema
x0 = [0; 5; 13; 17; 20];  % Condições iniciais do sistema

[y, t, x] = initial(t2, x0, t);  % Calcula a resposta do sistema com as condições iniciais

n_states = size(x, 2)  % Número de estados do sistema


figure
plot(t, x(:, 1), 'k', t, x(:, 2), 'c', t, x(:, 3), 'm', t, x(:, 4), 'y', t, x(:, 5), 'b');
legend('Estado 1', 'Estado 2', 'Estado 3', 'Estado 4', 'Estado 5');


%discreta
gz = c2d(sys, ta, 'zoh'); % Sistema Discreto em malha fechada
[az, bz, cz, dz, tz] = ssdata(gz); %Separa Variaveis

aez = [az zeros(2,1); -cz 1]; %Matriz A Expandida
bez = [bz; 0];  %Matriz B Expandida

z1 = exp(p1*ta); %Polos para tempo discreto
z2 = exp(p2*ta);
z3 = exp(p3*ta);

kez = acker(aez, bez, [z1 z2 z3]); %Ganhos em tempo Discreto
kz = kez(1:2); %Ganho da retroação
jz = -kez(3); %Ganho da Expansão

ys = ss([A-B*k1 B*jj; -C 0], 20*[0; 0; 1], [C 0]); %Matriz com Ganhos
us = ss([A-B*k1 B*jj; -C 0], [0; 0; 1], [-k1 jj]);

#expansão de polos explicita
yp = ss([az-bz*kz bz*jz; -cz 1], [0; 0; 1], [cz 0], 0, ta);
up = ss([az-bz*kz bz*jz; -cz 1], [0; 0; 1], [-kz jz], 0, ta);

#expansão de polos implicita
yc = ss([az-bz*kz bz*jz; -cz 1], [bz*jz; 1], [cz 0], 0, ta);
uc = ss([az-bz*kz bz*jz; -cz 1], [bz*jz; 1], [-kz jz], jz, ta);

%Tuste
yt = ss([az-bz*kz bz*jz; -cz 1], [bz*jz/2; 1], [cz 0], 0, ta);
ut = ss([az-bz*kz bz*jz; -cz 1], [bz*jz/2; 1], [-kz jz], jz/2, ta);

Lz = acker(az', (cz*az)', 1/3*[z1 z2]); %Ganho Ackermam // Polo desacelerado
a= bz*jz


Af = [az-bz*kz bz*jz bz*kz; -cz 1 0 0; [0;0] [0;0] [0;0] az-Lz'*cz]; % Matriz A - Retroação de estados observados com uma expansão de polos
Bf = [0; 0; 1; 0; 0]; % Matriz B - Retroação de estados observados com uma expansão de polos
Cf = [cz 0 0 0]; % Matriz C - Retroação de estados observados com uma expansão de polos

%saida
ypf= ss(Af,Bf*20,Cf,0,ta); %explicita
ycf = ss(Af,[bz*jz;1;0;0]*20,Cf,0,ta); %implicita
ytf = ss(Af,[bz*jz/2; 1; 0; 0]*20, Cf, 0, ta); %tustin

%entrada
upf = ss(Af, Bf, [-kz jz Lz], 0, ta); %explicita
ucf = ss(Af,[bz*jz;1;0;0], [-kz jz Lz], jz, ta); %implicita
utf = ss(Af,[bz*jz/2;1;0;0], [-kz jz Lz], jz/2, ta); %tustin
Lz
figure;
subplot(2, 1, 1);
step(ys, 2);
hold on;
step(ypf, 'r');
hold on;
step(ycf, 'g');
hold on;
step(ytf, 'k');
ylabel('Saida');

subplot(2, 1, 2);
step(us, 2);
hold on;
step(upf, 'r');
hold on;
step(ucf, 'g');
hold on;
step(utf, 'k');
ylabel('Entrada');
