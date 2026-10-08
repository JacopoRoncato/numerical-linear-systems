% Progetto 1 - Fondamenti di Calcolo Numerico (2026)
% Versione corretta per riproducibilita'' del codice pubblicato.
% Richiede nella stessa cartella: fwsub.m, bksub.m, jacobi.m,
% graddyn.m, gradprec.m, gc.m.
% Le implementazioni delle funzioni esterne non sono modificate qui.

clear
clc
close all

%% DATI DEL PROBLEMA
Vin = 15;
Rin = 60;
R=ones(6,1);
G=ones(6,1);

for i=1:6
    R(i) = Rin / (2*i);   % R1,...,R6
    G(i) = 1 / R(i);      % G1,...,G6
end
Gin = 1 / Rin;

G1 = G(1);
G2 = G(2);
G3 = G(3);
G4 = G(4);
G5 = G(5);
G6 = G(6);

%% PUNTO 2.1 - MATRICE A E TERMINE NOTO b

A= [ Gin+G1+G2,   -G1,         0,      -G2,         0,        0,        0,        0,        0;
         -G1,   2*G1+G2,      -G1,        0,      -G2,        0,        0,        0,        0;
           0,      -G1,      G1+G2,       0,        0,      -G2,        0,        0,        0;
         -G2,        0,         0,    G2+G3+G4,   -G3,        0,      -G4,        0,        0;
           0,      -G2,         0,      -G3,   G2+2*G3+G4,  -G3,        0,      -G4,        0;
           0,        0,       -G2,        0,      -G3,    G2+G3+G4,     0,        0,      -G4;
           0,        0,         0,      -G4,        0,        0,    G4+G5+G6,   -G5,        0;
           0,        0,         0,        0,      -G4,        0,      -G5,   G4+2*G5+G6,  -G5;
           0,        0,         0,        0,        0,      -G4,        0,      -G5,   G4+G5+G6 ];

b = [Gin*Vin; 0; 0; 0; 0; 0; 0; 0; 0];

disp('Matrice A =');
disp(A);

disp('Termine noto b =');
disp(b);

%% PUNTO 2.2 - VERIFICA ESISTENZA E UNICITA'' DELLA LU SENZA PIVOTING
% Condizione: tutti i minori principali devono essere non nulli
format short e

delta = zeros(9,1);
verificaLU = 0;

for k = 1:9
    delta(k) = det(A(1:k,1:k));
    if delta(k) == 0
        verificaLU = 1;
    end
end

disp('Determinanti di matrici minori =')
disp(delta)

disp('Esiste ed e'' unica la fattorizzazione LU con Lii=1 ?')
if verificaLU==1 disp('No')
else disp('Sì')
end

%% PUNTO 2.3 - FATTORIZZAZIONE LU CON MATLAB
[L,U,P] = lu(A);

detU = det(U);
disp('det(U) =')
disp(detU)

if norm(P - eye(9))==0 
    disp('Non è stato eseguito pivoting')
else disp('è stato eseguito pivoting')
end

%% PUNTO 2.4 - RISOLUZIONE CON fwsub E bksub
format long e

y  = fwsub(L, P*b);
xc = bksub(U, y);

disp('norma euclidea di xc =')
disp(norm(xc))

%% PUNTO 2.5 - RISOLUZIONE CON BACKSLASH MATLAB
xm = A\b;

disp('norma euclidea di xm =')
disp(norm(xm))

%% PUNTO 2.6
format short e

disp('norma infinito di (xm - xc) =')
disp(norm(xm - xc, inf))

%% PARTE 3
%%PUNTO 3.1
v=eig(A); dp=0;
for i=1:9 
    if (v(i)<=0) dp=1;
    end
end

if isequal(A, A') && all(v > 0)
    disp('A è simmetrica e definita positiva')
end

%%Punto 3.2
format long e

x0=zeros(9,1); toll=1e-11; nmax=1000;
[xJ,kJ]=jacobi(A,b,x0,toll,nmax);
disp('Iterazioni con metodo di Jacobi='); disp(kJ);

%%Punto 3.3
format short e

true_rel_err=norm(xm-xJ)/norm(xm);

%%Punto 3.4
D= diag(diag(A));
Bj= eye(size(A))-D\A;
p_J= max(abs(eig(Bj)));

confronto = p_J^kJ;
disp('err_relativo='); disp(true_rel_err);
disp('raggio spettrale^kJ'); disp(confronto);
%Otteniamo p_J^kJ= 1.4294e-11 e true_rel_err = 9.7673e-12
% Poiche' p_J = rho(B_J) < 1, il metodo di Jacobi converge.
% Inoltre, dalla teoria dei metodi iterativi stazionari, l'errore
% asintoticamente decresce con legge proporzionale a p_J^kJ.
% Pertanto il valore P_J^kJ fornisce una stima dell'ordine di
% grandezza dell'errore dopo kJ iterazioni.
% Il fatto che l'errore relativo calcolato risulti minore di P_J^kJ 
% conferma la coerenza dei risultati numerici con la teoria.

%%Punto 3.5
 [xG, iterG, errG] = graddyn (A, b, x0, nmax, toll);
 true_rel_errG= norm(xm-xG)/norm(xm);
 disp('Sono state necessarie, con il metodo del gradiente,');
 disp(iterG); disp('iterazioni e con un errore relativo pari a'); 
 disp(true_rel_errG);

 lambda_max=max(abs(eig(A))); lamda_min= min(abs(eig(A)));
 alpha_opt= 2/(lambda_max + lamda_min);
 Br=eye(size(A))-alpha_opt*A;
 p_optR=(cond(A)-1)/(cond(A)+1); %va bene come formula siccome abbiamo 
 %verificato che A è sdp

 confronto= p_optR^iterG;

 %kG             = 348
 %true_rel_err_G = 1.388966599318782e-11
 %rhoR           = 9.319942974997011e-01
 %confronto      = 2.268886497603655e-11


%//

%%Punto 3.6
P=diag(diag(A));
[xGP,iterGP,errGP] = gradprec(A,b,P,x0,nmax,toll);
true_rel_errGP= norm(xm-xGP)/norm(xm);

disp('Sono state necessarie, con il metodo del gradiente precondizionato,');
 disp(iterGP); disp('iterazioni e con un errore relativo pari a'); 
 disp(true_rel_errGP);

%%Punto 3.7
[xGC, iterGC, errGC] = gc (A, b, x0, nmax, toll);
true_rel_errGC= norm(xm-xGC)/norm(xm);

disp('Sono state necessarie, con il metodo del gradiente coniugato,');
 disp(iterGC); disp('iterazioni e con un errore relativo pari a'); 
 disp(true_rel_errGC);

%% Punto 3.8 - CONFRONTO DELLA CONVERGENZA
figure;
semilogy(0:iterG, errG, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Gradiente');
hold on;
semilogy(0:iterGP, errGP, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Gradiente Precondizionato');
semilogy(0:iterGC, errGC, 'b-.', 'LineWidth', 1.5, 'DisplayName', 'Gradiente Coniugato');
xlabel('Iterazioni (k)');
ylabel('||r^{(k)}|| / ||b||');
title('Confronto convergenza: Gradiente vs GP vs GC');
legend('show', 'Location', 'best');
grid on;
hold off;

%% Punto 3.9 - CONDIZIONAMENTO
K_A = cond(A);
K_PA = cond(P\A);
fprintf('Condizionamento di A: %.6e\n', K_A);
fprintf('Condizionamento di P^{-1}A: %.6e\n', K_PA);
