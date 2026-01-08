
function [E_max,i_max,j_max, k_max,E_rms_max,i_rms,j_rms, E_rms] = fiterrs_pos(k_red_test,Qjj_test,H_hat)
% E_max:          maximal difference between Qjj and H_hat at any panel any frequency
% E_rms_max:      maximal rms of error over all frequencies at any panel
% E_rms:          maximal rms of error over all frequencies and over all panels

[nj,nq,nk]     = size(Qjj_test);
E_hat = Qjj_test-H_hat;


E_rms = rms(E_hat,"all");

E_rms_jj = zeros(nj,nq);
for i = 1:nj
    for j = 1:nq
        E_rms_jj(i,j) = rms(E_hat(i,j,:));
    end
end
[E_rms_max,I_rms] = max(E_rms_jj,[],'all');
[i_rms,j_rms] = ind2sub(size(E_rms_jj),I_rms);
% disp(['Maximal rms over k error at entry ',num2str(i_rms),',',num2str(j_rms),'  :    ',num2str(E_rms_max*100),' %'])



% E_max
[E_max,I_max] = max(abs(E_hat),[],'all');
[i_max,j_max, k_max] = ind2sub(size(E_hat),I_max);
% disp(['Maximal error at k=',num2str(k_red_test(k_max)),' in entry  ',num2str(i_max),',',num2str(j_max),'  :    ',num2str(E_max*100),' %'])



end