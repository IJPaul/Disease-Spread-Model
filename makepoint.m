function P = makepoint(x,y,t)
% P is an infected person at a point with 
% x-coordinate P.x 
% y-coordinate P.y
% and t units left until being in the recovery state
P = struct('x', x,'y', y, 't',t);
end
