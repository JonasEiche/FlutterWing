function [PHI,OMEGA] = build_PHIgf(Kgg,Mgg,n_f)
% build_PHIgf       Structural modal analysis for normalized modeshapes and
%                   eigenfrequencies
%
% INPUT
%   Kgg             Global FEM DoF Stiffness Matrix
%   Mgg             Global FEM DoF Mass Matrix
%   nf              Number of evaluated Modes
%
% OUTPUT
%   PHI             [ng x nf] Modeshapes as column vectors from lowest
%                   eigenfrequency on the left to highest on the right
%   OMEGA           Vector of eigenfrequencies (rad) from low to high
%
% Jonas * July 2025
% _________________________________________________________________________


% Dirichlet Boundary Condition on the root node
Mgg_ = Mgg(4:end,4:end);
Kgg_ = Kgg(4:end,4:end);

if ~(n_f<=size(Mgg_,1))     % More Modes requested nf than DoFs after fixing left root
    warning(['Max number of modes is: ', num2str(size(Mgg_,1))])
    n_f = size(Mgg_,1);
end
% Dynamic Eigenmodes
try chol(Mgg_);
    % disp('Mass Matrix Mgg_ is symmetric positive definite.')
catch
    warning('Mass Matrix Mgg_ is not symmetric positive definite')
end
try chol(Kgg_);
    % disp('Stiffness Matrix Kgg_ is symmetric positive definite.')
catch
    warning('Stiffness Matrix Kgg_ is not symmetric positive definite')
end
[PHI,OMEGA_sq] = eigs(Kgg_,Mgg_,n_f,'smallestabs');

% [V,D] = eigs(A,B) returns V as a matrix whose columns are the generalized right eigenvectors that satisfy A*V = B*V*D. The 2-norm of each eigenvector is not necessarily 1.
% If B is symmetric positive definite, then the eigenvectors in V are normalized so that the B-norm of each is 1. If A is also symmetric, then the eigenvectors are B-orthonormal.


% sort modes according to frequency
[OMEGA_sq,ind] = sort(diag(OMEGA_sq));
PHI  = PHI(:,ind);
OMEGA = sqrt(OMEGA_sq);

% Normalization by setting generalized mass to unity is MATLAB default, hence is not necessary if Mgg is symmetric and positive definite, then PHI'*Mgg*PHI = eye is the standard normalization
% norm_gen_mass = diag(PHI' * Mgg_ * PHI);
% PHI = PHI ./ transpose(norm_gen_mass);

% % Normalize by setting the largest component of every modeshape vector to unity
idx_maxabs = max(PHI(:,:)) > -1*min(PHI(:,:));
norm_f = zeros(1,n_f);
norm_f(idx_maxabs)  = max(PHI(:,idx_maxabs));
norm_f(~idx_maxabs) = min(PHI(:,~idx_maxabs));
PHI = PHI./repmat(norm_f,length(PHI),1);

% Dirichlet Boundary Condition on the root node
PHI = [zeros(3,n_f);
       PHI];

end