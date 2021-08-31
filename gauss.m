clear
close all;

%normal matrix
%A=[1,2,3;4,5,6;7,8,9];
%A=A+diag([10,10,10]);
%x_e = [1,1,1]';

%big matrix
%A=[1,2,3,4,5,6;7,8,9,1,2,3;1,2,3,4,5,6;2,3,4,5,6,7;1,2,5,6,3,7;2,4,7,4,1,8];
%A = A +diag([10,10,10,10,10,10]);
%x_e = [1,1,1,1,1,1]';

%small matrix
%A=[1,2;4,5];
%A=A+diag([10,10]);
%x_e = [1,1]';

%diag identity matrix
%A=[1,2,3;4,1,6;7,8,1];
%x_e = [1,1,1]';

%diag 2 matrix
%A=[2,2,3;4,2,6;7,8,2];
%x_e = [1,1,1]';

%diag 3 matrix
%A=[3,2,3;4,3,6;7,8,3];
%x_e = [1,1,1]';

%diag 5 matrix
%A=[5,2,3;4,5,6;7,8,5];
%x_e = [1,1,1]';

%diag 11 matrix
%A=[1,2,3;4,1,6;7,8,1];
%A=A+diag([10,10,10]);
%x_e = [1,1,1]';

%symmetric matrix
%A=[1,2,4,5;2,3,5,6;4,5,6,7;5,6,7,8];
%A=A+diag([10,10,10,10]);
%x_e = [1,1,1,1]';

%symmetric matrix
%A=[1,2,4,5;2,3,5,6;4,5,6,7;5,6,7,8];
%x_e = [1,1,1,1]';

%random matrix
%rng('default');
%A=rand([0,10],[4,4]);
%A=A+diag([10,10,10,10]);
%x_e = [1,1,1,1]';

%sparse matrix
%A=[1,0,0,0,0;0,0,0,5,0;9,0,0,0,0;0,0,6,0,0;0,0,0,0,0];
%A=A+diag([10,10,10,10,10]);
%x_e = [1,1,1,1,1]';

%b = A*x_e;

nx = 80;
ny = nx;
sgrid = [nx,ny];
La = full(delsq(numgrid('S',nx+2)));

M = 40; %number of constraints
ix = sub2ind( sgrid, randi( nx, M, 1 ), randi( ny, M, 1 ));
%checker = 1;
% while checker == 1
%     update = 0;
%     for a=1:M
%         for b = a:M
%             if ix(a,:) == ix(b,:)
%                 update = 1;
%                 ix = sub2ind( sgrid, randi( nx, M, 1 ), randi( ny, M, 1 ));
%                 break
%             end
%         end
%         if update == 1
%             break
%         end
%     end
%     if update == 0
%         checker = 0;
%     end
% end
% i1 = sub2ind( sgrid, nx/4, ny/4 );
% i2 = sub2ind( sgrid, nx/4, ny/2 );
% i3 = sub2ind( sgrid, nx/2, ny/2 );
%J = zeros( M, size(La,1) );

 J = sparse( 1:M, ix, ones(M,1), M, size(La,1) );

 N = M;  % should rename M above to be N to avoid this conflict
%J(1,i1) = 1;
%J(2,i2) = 1;
%J(3,i3) = 1;
M = eye(size(La));
% h = 0.01;  % does larger h make convergence slower？ no?      
 %q = 1e5;       
 %A = J*((M+h^2*q*La)\J');
%sc = 1e-8;
sc = 1e-20;
A = J*((sc*M+La)\J');
%b=[10;-10;10];
%b = [1;-1;2];
%b=[53;-118;211];


b = (rand(N,1)-0.5) *4;
x_e = A\b;

x = zeros(N,1);
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
eps=1e-40;
maxiter = 500;
deltaXNorm = zeros(1,maxiter);
errorNorm = zeros(1,maxiter);
error = zeros(1,maxiter);
errorUp = zeros(1,maxiter);
errorDown = zeros(1,maxiter);
errorSum = zeros(1,maxiter);
y = x;

xLo = zeros(N,1);
xHi = zeros(N,1);
%xHi(1:3:end) = Inf;
xLo(1:end)  = -2;
xHi(1:end) = 2;
%x_e = min(max(x_e,xLo),xHi)';
u_t = (sc*M+La)\(J');
for k = 1:maxiter
    F = [];
    T = [];
    xbound = zeros(N,1);
    fInd = 1;
    tInd = 1;
    for i=1:N
        x(i) = (b(i)-L(i,:)*x-U(i,:)*x)/A(i,i);
        %x(i) = x(i) + (b(i) - A(i,:)*x)/A(i,i);
        if xLo(i) >= x(i)
            T(tInd) = i;
            tInd = tInd + 1;
            xbound(i) = xLo(i); 
        elseif xHi(i) <= x(i)
            T(tInd) = i;
            tInd = tInd + 1;
            xbound(i) = xHi(i);
        else
            F(fInd) = i;
            fInd = fInd +1;
        end
        x(i) = min( max(xLo(i),x(i)), xHi(i));
    end
%     if norm(x-y) < eps
%         break;
%     end
    deltaXNorm(k) = norm(x-y);
    errorNorm(k) = norm(x-x_e);  % not correct with bounds
    
    fLength = length(F);
    tLength = length(T);
    sum = 0;
    for fLoop=1:fLength
        for fsLoop=1:fLength
        sum = sum + A(F(fLoop),F(fsLoop))*x(F(fLoop));
        end
        sum = sum - b(F(fLoop));
    end
    for fLoop = 1:fLength
        for tLoop = 1:tLength
            sum = sum + A(F(fLoop),T(tLoop))*xbound(T(tLoop));
        end
    end
    %error(k) = norm(A*x-b); % F=[1,2,4]; A(F,F) * x(F) - b(F) + A(F,T) * xbound(T)
    error(k) = norm(sum);
    
    matrixUpSum = J'*x_e;
%     for u_oneLoop = 1:N
%         if x(u_oneLoop) == xLo(u_oneLoop) 
%             matrixSum = matrixSum - J(u_oneLoop,:)'*xLo(u_oneLoop);
%         elseif x(u_oneLoop) == xHi(u_oneLoop)
%             matrixSum = matrixSum - J(u_oneLoop,:)'*xHi(u_oneLoop);
%         else
%             matrixSum = matrixSum - J(u_oneLoop,:)'*x(u_oneLoop);
%         end
%     end
    for u_tLoop = 1:N
        matrixUpSum = matrixUpSum - J(u_tLoop,:)'*x(u_tLoop);
    end
    %matrixUpSum = matrixUpSum - J'*x;
    errorUp(k) = norm(matrixUpSum);
    
    
%     matrixDownSum = 0;
%     for u_downLoop = 1:N
%         if x(u_downLoop) ~= xLo(u_downLoop) && x(u_downLoop) ~= xHi(u_downLoop)
%             matrixDownSum = matrixDownSum + J(u_downLoop)*x(u_downLoop) - 1;
%         end
%     end
    
    %matrixDownSum = 0;
   % downS = u_t*x;
%     for u_downLoop = 1:N
%         if x(u_downLoop) ~= xLo(u_downLoop) && x(u_downLoop) ~= xHi(u_downLoop)
%             matrixDownSum = matrixDownSum + J(u_downLoop)*downS - 1;
%         end
%     end
    matrixDownSum =  J*u_t*x - b;
    %matrixDownSum = (inv(sc*M+La)*(J'));
    errorDown(k) = norm(matrixDownSum);
    
    errorSum(k) = errorUp(k)+ errorDown(k);
    
    y = x;
end
x;
norm(x-x_e);

%u_o = b/J*inv(M+h^2*q*La);
%u_t = inv(M+h^2*q*La)*(u_o + J'*x);

figure
semilogy(deltaXNorm);
xlabel('number of iteration');
ylabel('change in x');
% figure
% semilogy(errorNorm);
% xlabel('number of iteration');
% ylabel('error value');

figure
semilogy(error);
xlabel('number of iteration');
ylabel('error in x');

figure
semilogy(errorUp(1:50));
xlabel('number of iteration');
ylabel('upper part of matrix error in x');

figure
semilogy(errorDown(1:50));
xlabel('number of iteration');
ylabel('lower part of matrix error in x');

figure
semilogy(errorSum(1:50));
xlabel('number of iteration');
ylabel('whole matrix error in x');

figure
%xcontraint = (M+h^2*q*La)\(J'*x);
xcontraint = (sc*M+La)\(J'*x);
imagesc( reshape( xcontraint, sgrid ));
colorbar;
hold
for i=1:N
    if x(i) == xLo(i) || x(i) == xHi(i)
        jInd = find(J(i,:)==1);
        [row,col] = ind2sub(sgrid,jInd);
        plot(col,row, 'ro', 'MarkerSize', 10);
    end
end

figure
[X,Y] = meshgrid(1:sgrid(1),1:sgrid(2));
Z = reshape(xcontraint,sgrid);
C = Z;
surf(X,Y,Z,C);

%LPL = (M+h^2*q*La)\(J'*x_e);
%imagesc( reshape( LPL, sgrid ));
%colorbar;

%figure
%[X,Y] = meshgrid(1:sgrid(1),1:sgrid(2));
%Z = reshape(LPL,sgrid);
%C = Z;
%surf(X,Y,Z,C);
