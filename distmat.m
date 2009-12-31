function [ D ] = distmat( v )
%DISTMAT Finds the distances between every two elements in the vector

D = [];
for i = 1 : (length(v)-1)
    for j = i+1 : length(v)
        d = abs(v(i) - v(j));
        D(i,j) = d;
        D(j,i) = d;
    end
end

end

