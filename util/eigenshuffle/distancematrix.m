function d = distancematrix(vec1,vec2)

% simple interpoint distance matrix
[vec1,vec2] = ndgrid(vec1,vec2);
d = abs(vec1 - vec2);

end