function gramDet = gramDet(A,varargin)
%GRAMDET  Gram determinant of a set of column vectors.
%   D = GRAMDET(A) returns det(G) where G = A'*A.
%
%   D = GRAMDET(A,'normalize',true) first scales each column of A to
%   unit-length, making the determinant invariant to column scaling.
%
%   The determinant is computed via the singular values of A for
%   numerical robustness (det(A'*A) = prod(s.^2)).
%
%   Examples
%   --------
%   % 4 random 10-D vectors, scale-invariant volume
%   A = randn(10,4);
%   d = gramDet(A,'normalize',true)
%
%   % Compare with the raw, scale-dependent determinant
%   detRaw = gramDet(A)
%
%   See also VECNORM,  SVD,  DET.

% Jonas * 08.2025 ----------------------------------------------

%--- Parse optional arguments -------------------------------------------
p = inputParser;
p.addParameter('normalize',false,@islogical);
p.parse(varargin{:});
doNorm = p.Results.normalize;

%--- Normalise columns if requested -------------------------------------
if doNorm
    colNorms = vecnorm(A,2,1);           % fast column norms :contentReference[oaicite:0]{index=0}
    % handle zero columns (volume must be zero)
    if any(colNorms < 1e-8)
        gramDet = 0;
        return;
    end
    A        = A ./ colNorms;
end

%--- If more vectors than ambient dimension: determinant is zero
[m,n] = size(A);
if m < n
    gramDet = 0;
    return;
end

%--- Compute determinant through singular values ------------------------
%   For an m×n matrix (m>=n assumed), det(A'*A) = prod(σ_i^2).
s  = svd(A,'econ');                      % economical SVD :contentReference[oaicite:1]{index=1}
gramDet = prod(s.^2);

%--- Force real output (should already be real within round-off) --------
gramDet = real(gramDet);
end
