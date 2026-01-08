function [reorder,e1,phi1] = eigenshuffle(e0,phi0,e1,phi1)
% taken from John D'Errico eigenshuffle.m function

dist = (1-abs(phi0'*phi1)).*sqrt( ...
    distancematrix(real(e0),real(e1)).^2+ ...
    distancematrix(imag(e0),imag(e1)).^2);
  
  reorder = munkres(dist);
  
  phi1 = phi1(:,reorder);
  e1   = e1(reorder);
  
  S = squeeze(real(sum(phi0.*phi1,1))) < 0;
  phi1(:,S) = -phi1(:,S);
  
end