% verify derivative formula of Eig and SVD
fnormalize = @(x) x ./ sqrt(real(sum(conj(x).*x, 1)));

rand('state', 1312);
n = 5;

%%%%%%%%%%%%%%%%%%%%%
% Real symmetric case
de = 1e-4;
A = rand(n);
A = (A + A')/2;
dA = de * rand(n);
%dA = (dA + dA')/2;

[U, D] = eigh(A);
[Up, Dp] = eigh(A + dA);

k = 3;

dform_lambda = U(:,k)' * dA * U(:,k)
delta_lambda = Dp(k) - D(k)

dform_vec = pinv(D(k) * eye(n) - A) * dA * U(:,k)
delta_vec = Up(:,k) - U(:,k)

%%%%%%%%%%%%%%%%%%%%%
% general case, diagC=0
de = 1e-4;
A = rand(n);
dA = de * rand(n);

[U, D] = eigh(A);
[Up, Dp] = eigh(A + dA);

dform_lambda = diag(inv(U) * dA * U)
delta_lambda = Dp - D

F = 1 ./ ((D - D.').' + Inf * eye(n));

dform_vec = U * (F .* (U \ dA * U));
dform_vec = fnormalize(U + dform_vec) - U;
delta_vec = Up - U;

dform_vec(:,1)
delta_vec(:,1)