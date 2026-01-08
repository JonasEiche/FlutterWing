function [poles,A0,A1,A2,Ap,D,E,R] = rogersRFA_magW(k_red, Qjj, np)
%   Roger’s rational-function approximation with weighting by entry magnitude.
%   A0 = Qjj(k_red=0) !
%   A2 = 0 ! 
%
%   INPUTS
%     k_red  1×nk        vector of reduced frequencies, *first entry must be 0*
%     Qjj    ny×nu×nk    array of AIC matrices evaluated at k_red
%     np     scalar      desired number of lag poles 
%
%   OUTPUTS
%     poles  1×np        positive pole magnitudes p_r  (R uses -p_r)
%     A0     ny×nu       quasi-steady term  (k = 0)
%     A1     ny×nu       proportional to i·k
%     A2     ny×nu       proportional to −k²  (unused but returned)
%     AP     ny×nu×np    lag-state gain matrices A_r
%     D,E,R  —           state-space realisation of the input side
%            x_dot  = R x + E u,   y = D x + A0 u + A1 u_dot
% 
%   ALGORITHM
%   ---------
%   Roger’s Rational Function Approximation writes each element Q̂(k) as
%
%            Q̂(k) = A0 + ik A1  − k^2 A2
%                    + \sum_{r=1}^{np} A_r ik/(ik + p_r)      (1)
% 
%                 = A0 + ik A1  − k^2 A2
%                    + \sum_{r=1}^{np} (A_r ik p_r + A_r k^2)/(k^2 + p_r^2).
%
%   After subtracting A0 the unknown real matrices [A1 , A_r] are obtained
%   from a *real* least-squares fit of (1) to the supplied data.
%   Because A_r appear in both real and imaginary parts we build, for every
%   k≠0, the 2×(np+1) block.
%
%           ┌────────────┬────────────────────────────────────┐
%           │ 0          │      k²/(k²+p_r²)  (r = 1…np)      │   ← Re
%       A = │────────────┼────────────────────────────────────│
%           │ k          │  k p_r/(k²+p_r²)  (r = 1…np)       │   ← Im
%           └────────────┴────────────────────────────────────┘          
%
%   Stacking the real & imaginary parts of Q̂(k) on the RHS gives an
%   over-determined real system A·x = b that is solved by x = A\b.
%
%
% Jonas * July 2025
% _________________________________________________________________________

[ny,nu,nk] = size(Qjj);
assert(nk==length(k_red),'number of reduced frequencies must match AIC Matrices')
poles = max(k_red)./(np:-1:1);      % [kmax/6  kmax/5  kmax/4  kmax/3  kmax/2  kmax/1]


% A0=Q(0),A1
assert(k_red(1)==0,'first reduced frequency must be zero. Quasisteady VLM contribution')
A0 = Qjj(:,:,1);
A2 = zeros(ny,nu);
nk_=nk-1;
Qjj_ = Qjj(:,:,2:end)-A0;
k_red_ = k_red(2:end);

A_lin = zeros(2*nk_,np+1); 
for ii=1:nk_
  A_lin(2*ii-[1 0],:) = [0,            k_red_(ii)^2./( k_red_(ii)^2 + poles.^2 );       % real
                         k_red_(ii),   k_red_(ii)*poles./( k_red_(ii)^2 + poles.^2 )];  % imag
end

W(:,:,(1:2:2*nk)+1) = 1./max(abs(Qjj),1e-10); % real weight = imag weight
W(:,:,(1:2:2*nk))   = 1./max(abs(Qjj),1e-10); % imag weight = real weight
W = W(:,:,3:end);


B_lin = zeros(ny,nu,np+1);
for jj=1:nu % columns/input
    for ii=1:ny % rows/output
        b_lin = reshape([reshape(real(Qjj_(ii,jj,:)),1,nk_,1);reshape(imag(Qjj_(ii,jj,:)),1,nk_,1);],2*nk_,1,1);
        w = diag(reshape(W(ii,jj,:),2*nk_,1,1));
        B_lin(ii,jj,:) = reshape(  (w*A_lin)\(w*b_lin),   1,1,[]);
    end
end

A1 = B_lin(:,:,1);
Ap = B_lin(:,:,2:end);

% Input Realization with nu*np lag states
% D      = reshape(B_lin(:,:,2:end),ny,nu*np);
D      = reshape(Ap,ny,nu*np);
% % E      = repmat(speye(nu),np,1);
% % R      = diag(sparse(reshape(repmat(-poles, nu,1),1,nu*np)));
E      = repmat(eye(nu),np,1);
R      = diag(reshape(repmat(-poles, nu,1),1,nu*np));


% % Output Realization with ny*np lag states
% D_r = zeros(ny*np,nu);
% for i=1:np
%     D_r((i-1)*ny+1:i*ny,:) = Ap(:,:,i);
% end
% E_r = repmat(eye(ny),1,np);
% R_r = diag(reshape(repmat(-poles, ny,1),1,ny*np));



end