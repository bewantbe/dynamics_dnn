% verify derivative formula of Eig and SVD
fnormalize = @(x) x ./ sqrt(real(sum(conj(x).*x, 1)));

rand('state', 1312);
n = 5;
m = 3;
de = 1e-4;

if 0
%%%%%%%%%%%%%%%%%%%%%
disp('Real symmetric case:');
A = rand(n);
A = (A + A')/2;
dA = de * rand(n);
%dA = (dA + dA')/2;

[U, D] = eigh(A);
[Up, Dp] = eigh(A + dA);

k = 3;

assert(sum(D == D(k)) == 1);  % simple root

dform_lambda = U(:,k)' * dA * U(:,k)
delta_lambda = Dp(k) - D(k)

dform_vec = pinv(D(k) * eye(n) - A) * dA * U(:,k)
delta_vec = Up(:,k) - U(:,k)

%%%%%%%%%%%%%%%%%%%%%
disp('General case, u0*u=1:')
A = rand(n);
dA = de * rand(n);

[U, D] = eigh(A);
[Up, Dp] = eigh(A + dA);

k = 3;

assert(sum(D == D(k)) == 1);  % simple root

[V, Dv] = eigh(A');

dform_lambda = V(:,k)' * dA * U(:,k) / (V(:,k)' * U(:,k))
delta_lambda = Dp(k) - D(k)

dform_vec = pinv(D(k) * eye(n) - A) * (eye(n) - (U(:,k) * V(:,k)') / (V(:,k)' * U(:,k))) * dA * U(:,k)
delta_vec = Up(:,k) - U(:,k)
%norm(U(:,k)+dform_vec)

end
%%%%%%%%%%%%%%%%%%%%%
disp('General case, diagC=0:')
rand('state', 1);
A = rand(n);
dA = de * rand(n);

[U, D] = eigh(A);
[Up, Dp] = eigh(A + dA);

assert(length(unique(D)) == length(D));  % no repeat roots

dform_lambda = diag(U \ dA * U)
delta_lambda = Dp - D

F = 1 ./ ((D - D.').' + Inf * eye(n));

k = 3;
dform_vec = U * (F .* (U \ dA * U));
%norm(U(:,k) + dform_vec(:,k))
% renormalize to match MATLAB
dform_vec = fnormalize(U + dform_vec) - U
delta_vec = Up - U

%disp('eigenvectors')
%dform_vec(:,k)
%delta_vec(:,k)

return

%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Derivative of SVD:')
A = rand(m, n);
dA = de * rand(m, n);

[U, S, V] = svd(A);
[Up, Sp, Vp] = svd(A + dA);

dform_sv = diag(U.' * dA * V)
delta_sv = diag(Sp - S)

dP = U.' * dA * V;
dP1 = dP(:,1:m);
dP2 = dP(:,m+1:n);

S1 = S(:,1:m);

vS = diag(S);
F = 1 ./ ((vS.^2 - (vS.^2).').' + Inf * eye(m));

dform_u = U * (F .* (dP1 * S1 + S1 * dP1.'))
delta_u = Up - U

dform_D2 = S1 \ dP2;
dform_D1 = F .* (S1 * dP1 + dP1.' * S1);

dform_v2f = V * [dform_D1; dform_D2.']
delta_v = Vp - V

