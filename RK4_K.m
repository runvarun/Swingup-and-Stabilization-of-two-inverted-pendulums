%% RK 4 for K
function X = RK4_K(A,B,E,X0,T,R,Q)
   X=zeros(size(A,1),size(A,2),size(A,3));
   X(:,:,end)=X0;
   h=T/size(A,3)/2;
   for i=size(A,3):-1:2
        f1 = RHS(X(:, :, i), E(:, :, i), A(:, :, i), B(:, :, i), R, Q);
        f2 = RHS(X(:, :, i)-f1/2*h, E(:, :, i), A(:, :, i), B(:, :, i), R, Q);
        f3 = RHS(X(:, :, i)-f2/2*h, E(:, :, i), A(:, :, i), B(:, :, i), R, Q);
        f4 = RHS(X(:, :, i)-f3*h, E(:, :, i), A(:, :, i), B(:, :, i), R, Q);
        X(:, :, i-1) = X(:, :, i) - h*(f1/6+((f2+f3)/3)+f4/6);
   end

   function dX=RHS(X,E,A,B,R,Q)
        dX = -(E'^-1) * (A'*X*E + E'*X*A - E'*X*B*R^-1*B'*X*E + Q ) * E^-1;
   end
end