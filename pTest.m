function [P U V]  = pTest(X,Y,R,L)

dx = Y(2,1) - Y(1,1);

P = zeros(size(X));
U = zeros(size(X));
V = zeros(size(X));

%%%%%%%%%%%%%%The value of E might change%%%%%%
E = .5 * R;
%E = dx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = sqrt( (X-(L+R)).^2 + Y.^2 ) - R;

id = -E <= z & z <= E;
id2 = z < -E;

P(id) = -1/(2*R)*(1 - z(id)/E - sin(pi*z(id)/E)/pi);
P(id2) = -1/(R);

end