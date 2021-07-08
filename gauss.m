nx = 60;
ny = nx;
sgrid = [nx,ny];
La = full(delsq(numgrid('S',nx+2)));
i1 = sub2ind( sgrid, nx/4, ny/4 );
i2 = sub2ind( sgrid, nx/4, ny/2 );
i3 = sub2ind( sgrid, nx/2, ny/2 );
J = zeros(3,size(La,1));
J(1,i1) = 1;
J(2,i2) = 1;
J(3,i3) = 1;
M = eye(size(La));
h = 0.01;         
q = 1e5;       
A = J*((M+h^2*q*La)\J');
b=[1;-1;2];
x_e = A\b;

N = length(b);
x = zeros(N,1);
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
eps=1e-7;
deltaXNorm = zeros(1,10000);
errorNorm = zeros(1,10000);
y = x;

xLo = zeros(N,1);
xHi = zeros(N,1);
%xHi(1:3:end) = Inf;
xLo(1:end)  = -50;
xHi(1:end) = 50;

for k=1:10000
    for i=1:N
        x(i) = (b(i)-L(i,:)*x-U(i,:)*x)/A(i,i);
        x(i) = min(max(xLo(i),x(i)),xHi(i));
    end
    if norm(x-y) < eps
        break;
    end
        y = x;
end
x;
norm(x-x_e);



xcontraint = (M+h^2*q*La)\(J'*x);
imagesc( reshape( xcontraint, sgrid ));
colorbar;

figure
[X,Y] = meshgrid(1:sgrid(1),1:sgrid(2));
Z = reshape(xcontraint,sgrid);
C = Z;
surf(X,Y,Z,C);

%LPL = (M+h^2*q*La)\(J'*x_e);
%imagesc(reshape( LPL, sgrid ));
%colorbar;

%figure
%[X,Y] = meshgrid(1:sgrid(1),1:sgrid(2));
%Z = reshape(LPL,sgrid);
%C = Z;
%surf(X,Y,Z,C);
