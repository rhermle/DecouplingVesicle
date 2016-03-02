function [ F1 F2 ] = force(X, Y, mu, width, p0, g, R, dx, L)

graph = 0;

%E = dx;
E = .5 * R;
z = sqrt((X - (R + L)).^2 +Y.^2) - R;

id = (z >= -E) & (z <= E);

delta = zeros(size(X));
delta (id) = (1 + cos((pi * z(id) / E))) / (2 * E);

F1 = (delta/R) .* ((X - (R + L))) ./ sqrt((X - (R + L)).^2 + Y.^2 + eps);
F2 = (delta/R) .* (Y ./ (sqrt((X - (R + L)).^2 + Y.^2 + eps)));

% f1x = zeros(size(X));
% f1y = zeros(size(X));
% 
% f1x(id) = (1/(2*R*E)) * ((X(id) - (L+R)) ./ (sqrt( (X(id)-(L+R)).^2 + Y(id).^2 )) + (cos(pi*z(id)/E)) .*(X(id) - (L+R)) ./ (sqrt( (X(id)-(L+R)).^2 + Y(id).^2 + eps)));
% 
% f1y(id) = (1/(2*R*E)) * (Y(id) ./ (sqrt( (X(id)-(L+R)).^2 + Y(id).^2 )) + (cos(pi*z(id)/E)) .*Y(id) ./ (sqrt( (X(id)-(L+R)).^2 + Y(id).^2 + eps)));
% 
% divF = f1x + f1y;

if graph

    figure(12345);
    hold on;
    quiver(X,Y,F1,F2);
    title('F1 and F2');
    %quiver(x(1:quivRes:end,1:quivRes:end),Y(1:quivRes:end,1:quivRes:end),U(1:quivRes:end,1:quivRes:end, i),V(1:quivRes:end,1:quivRes:end, i));
    contour(X,Y,z,[0,0]);
    hold off;
    
end


end