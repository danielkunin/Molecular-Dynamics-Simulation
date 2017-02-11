function [ dist ] = avg_dist( pos, connect )
%   pos = 3 x n matrix of current positions (m)
%   connect = n x n array representing initial bond distance
%   ----------------------------------------------
%   dist = average dist between all pairs of bound atoms (m)

n = size(pos,2);
dist = 0;
num = 0;
for i = 1:n
    for j = i+1:n
        if connect(i,j)
            rij = norm(pos(:,i) - pos(:,j));
            dist = dist + rij;
            num = num + 1;
        end
    end
end
dist = dist / num;

end

