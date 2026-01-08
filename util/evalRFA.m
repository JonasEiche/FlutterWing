function [E_max,E_rms_max, E_rms, H] = evalRFA(k_red,Qjj,Q0jj,Q1jj,C,B,A)

[ny,nu,nk] = size(Qjj);

H = zeros(ny,nu,nk);
for i_sk = 1:nk
    sk = k_red(i_sk);
    H(:,:,i_sk) = C* ((1i*sk*eye(size(A,1))-A) \ B)*1i*sk + Q0jj + 1i*sk*Q1jj;
end
[E_max,E_rms_max, E_rms] = fiterrs(k_red,Qjj,H);

end