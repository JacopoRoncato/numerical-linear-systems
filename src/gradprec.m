function [x,iter,err] = gradprec(A,b,P,x0,nmax,toll)

% Metodo del gradiente precondizionato

% A: matrice del sistema lineare Ax = b
% b: termine noto
% P: precondizionatore (matrice) 
% x0: guess iniziale
% nmax: numero massimo di iterazioni
% toll: tolleranza sul residuo normalizzato

% x: soluzione ottenuta
% iter: numero di iterazioni effettuate
% err: vettore contenente il residuo normalizzato ad ogni iterazione

% Verifico che A matrice quadrata
[n,m] = size(A);
if(n~=m)
    error('matrice non quadrata')
end

% Inizializzo x a x0, numero di iterazioni, definisco residuo ed
% errore (residuo normalizzato)
x = x0;
iter = 0;
r = b - A*x;
errk = norm(r)/norm(b);
err = errk;

while(errk > toll && iter < nmax)
    iter = iter + 1;
    z = P\r;
    alpha = (z'*r) / (z'*A*z);
    x = x + alpha*z;
    r = r - alpha*A*z;
    errk = norm(r) / norm(b);
    err = [err errk];
end
if(iter == nmax)
    fprintf('il metodo non converge in %d iterazioni \n', iter);
end