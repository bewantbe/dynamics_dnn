% Reordered version of eig.
% Mainly for hermite matrix.

function [evecs, evals] = eigh(A)
  if nargout == 2
    [evecs, evals] = eig(A);
    evals = diag(evals);
    [evals, ord] = sort(evals, 'descend');
    evecs = evecs(:,ord);
    n = size(A,1);
    % flip the eigen vectors so that mostly positive.
    sign_v = sign(sum(evecs .* abs(evecs), 1));
    evecs = evecs .* sign_v;
  else
    evals = eig(A);
    evals = sort(evals, 'descend');
    evecs = evals;
  end

