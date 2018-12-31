% Reordered version of eig().
% Larger eigenvalues sort first.
% Flip eigenvectors so that they have more large positive elements.

function [evecs, evals] = eigh(A)
  if nargout == 1
    evecs = sort(eig(A), 'descend');
  else
    [evecs, evals] = eig(A);
    [evals, ord] = sort(diag(evals), 'descend');
    evecs = evecs(:,ord);
    % flip the eigen vectors so that mostly positive.
    sign_v = sign(sum(evecs .* abs(evecs), 1));
    evecs = evecs .* sign_v;
  end

