%% RK 4 for K
function P = RK4_L(A,C,E,P0,T,R,Q)
   P=zeros(size(A,1),size(A,2),size(A,3));
   P(:,:,1)=P0;
   h=0.01;
   for i=1:size(A,3)-1
        f1 = LHS(P(:, :, i), E(:, :, i), A(:, :, i), C(:, :, i)', R, Q);
        f2 = LHS(P(:, :, i)+f1/2*h, E(:, :, i), A(:, :, i), C(:, :, i)', R, Q);
        f3 = LHS(P(:, :, i)+f2/2*h, E(:, :, i), A(:, :, i), C(:, :, i)', R, Q);
        f4 = LHS(P(:, :, i)+f3*h, E(:, :, i), A(:, :, i), C(:, :, i)', R, Q);
        P(:, :, i+1) = P(:, :, i) + h*(f1/6+((f2+f3)/3)+f4/6);
   end

   function dP=LHS(P, E,A,C,R,Q)
        dP = (E'^-1) * (A'*P*E + E'*P*A - E'*P*C'*R^-1*C*P*E + Q ) * E^-1;
   end
end