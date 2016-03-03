function [ P U V X Y] = Stokes2DG(width, height, R, L, g, M, p0, mu, toGraph )

%2D Laplace
% p(0,y) = p0, p(1,y) = 100, p_y(x,0) = 0, p_y(x,1) = 0

% Ap=b

%Vectors p and A are indexed as:
% 1 
% 2
% 3
% 4 
% 5
% 6
% 7
% 8
% 9
% 10
% 11
% 12
% 13
% 14
% 15
% 16

%p:
% 1,1 
% 2,1
% 3,1
% 4,1
% 1,2 
% 2,2
% 3,2
% 4,2
% 1,3 
% 2,3
% 3,3
% 4,3
% 1,4  
% 2,4
% 3,4
% 4,4

%The following description is for M=4:
%A and p are indexed from 1 to 16 based on the subscripts for p above, in
%order to create 16 different equations.
%The second index for A is another 1 to 16 based on the same subscripts so
%each equation (row) for A matches up with the elements of p.

%p(i-1,j) - 4p_(i,j) + p(i+1,j) + p(i,j+1) + p(i,j-1) = 0


%M=30;

x = linspace(0,width,M);
y = linspace(-height/2, height/2 , M);

[X,Y] = meshgrid(x,y);

d = Y(2,1) - Y(1,1);

[F1 F2] = force(X, Y, mu, width, p0, g, R, d, L);

count = 1;
ind = zeros(M,M);

for i=1:M
    for j=1:M
        ind(i,j) = count;
        count = count + 1;
    end
end

%A = zeros(M*M,M*M);
A = sparse([],[],[],M*M,M*M,5*M*M);
b = zeros(M*M,1);

%interior points
for i=2:M-1
    for j=2:M-1
        
        index = ind(i,j);
        
        A(index, ind(i-1,j)) = 1;
        A(index, index) = -4;
        A(index, ind(i+1,j)) = 1;
        A(index, ind(i,j+1)) = 1;
        A(index, ind(i,j-1)) = 1;

        b(index) = d * (F1(j,i+1) - F1(j,i-1)) / (2) + d * (F2(j+1,i) - F2(j-1,i)) / (2);        
    end
end

%border conditions

for j = 1:M
    index1 = ind(1,j);
    indexM = ind(M,j);
    
    A(index1, index1) = 1;
    A(indexM, indexM) = 1;
    b(index1) = 0;
    b(indexM) = 0;
end

for i = 2:M-1
    index1 = ind(i,1);
    indexM = ind(i,M);
    
    A(index1, index1) = 1;
    b(index1) = 0;
    
    A(indexM, indexM) = 1;
    b(indexM) = 0;
end

A = sparse(A);
p = A\b;
%rcond(A) OR use cond

P = p(ind)';
    
if(toGraph)
    figure(1);
    surf(X,Y,P);
    title('Pressure');
    xlabel('-1 < x < 1');
    ylabel('-1 < y < 1');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Poisson Equation:
% u_xx + u_yy = p_x / mu
% u_x(0,y) = 0
% u_x(1,y) = 0
% u(x,0) = 0
% u(x,1) = 0

%Vectors u and A are indexed as:

% 1 
% 2
% 3
% 4 
% 5
% 6
% 7
% 8
% 9
% 10
% 11
% 12
% 13
% 14
% 15
% 16

%p:
% 1,1 
% 2,1
% 3,1
% 4,1
% 1,2 
% 2,2
% 3,2
% 4,2
% 1,3 
% 2,3
% 3,3
% 4,3
% 1,4  
% 2,4
% 3,4
% 4,4


%The following description is for M=4:
%A is first indexed from 1 to 16 based on the subscripts for u above, in
%order to create 16 different equations.
%The second index for A is another 1 to 16 based on the same subscripts so
%each equation (row) for A matches up with the elements of u.

%u(i-1,j) - 4u_(i,j) + u(i+1,j) + u(i,j+1) + u(i,j-1) = d*(p(i+1,j)-p(i-1,j))/2

%M=30;

%A = zeros(M*M,M*M);
A = sparse([],[],[],M*M,M*M,5*M*M);
b = zeros(M*M,1);

%interior points
for i=2:M-1
    for j=2:M-1
        
        index = ind(i,j);
        
        A(index, ind(i-1,j)) = 1;
        A(index, ind(i,j)) = -4;
        A(index, ind(i+1,j)) = 1;
        A(index, ind(i,j+1)) = 1;
        A(index, ind(i,j-1)) = 1;

        b(index) = d * (p(ind(i+1,j)) - p(ind(i-1,j)))/(2 * mu) - (d^2 / mu) * F1(j,i);      
    end
end

%border conditions

% u(x,0) = 0
% u(x,1) = 0
for i = 1:M
    index1 = ind(i,1);
    indexM = ind(i,M);
    
    A(index1, index1) = 1;
    A(indexM, indexM) = 1;
    b(index1) = 0;
    b(indexM) = 0;
end

for j = 2:M-1
    index1 = ind(1,j);
    indexM = ind(M,j);
    
    A(index1, index1) = 1;
    b(index1) = 0;
    
    A(indexM, indexM) = 1;
    b(indexM) = 0;
end

A = sparse(A);
u = A\b;

U = u(ind)';

if(toGraph)
    figure(2);
    surf(X,Y,U);
    title('U (Horizontal Velocity)');
    xlabel('-1 < x < 1');
    ylabel('-1 < y < 1');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Poisson Equation:
% v_xx + v_yy = p_y / mu
% v(0,y) = 0
% v(1,y) = 0
% v(x,0) = 0
% v(x,1) = 0

%M=30;

%A = zeros(M*M,M*M);
A = sparse([],[],[],M*M,M*M,5*M*M);
b = zeros(M*M,1);

%interior points
for i=2:M-1
    for j=2:M-1
        
        index = ind(i,j);
        
        A(index, ind(i-1,j)) = 1;
        A(index, ind(i,j)) = -4;
        A(index, ind(i+1,j)) = 1;
        A(index, ind(i,j+1)) = 1;
        A(index, ind(i,j-1)) = 1;

        b(index) = d * (p(ind(i,j+1)) - p(ind(i,j-1)))/(2 * mu) - (d^2 / mu) * F2(j,i);
    end
end

%border conditions

% v(x,0) = 0
% v(x,1) = 0
for i = 1:M
    index1 = ind(i,1);
    indexM = ind(i,M);
    
    A(index1, index1) = 1;
    A(indexM, indexM) = 1;
    b(index1) = 0;
    b(indexM) = 0;
end

% v(0,y) = 0
% v(1,y) = 0
for j = 1:M
    index1 = ind(1,j);
    indexM = ind(M,j);
    
    A(index1, index1) = 1;
    A(indexM, indexM) = 1;
    b(index1) = 0;
    b(indexM) = 0;
end

A = sparse(A);
v = A\b;

V = v(ind)';

if(toGraph)
    figure(3);
    surf(X,Y,V);
    title('V (Vertical Velocity)');
    xlabel('-1 < x < 1');
    ylabel('-1 < y < 1');
end

end

