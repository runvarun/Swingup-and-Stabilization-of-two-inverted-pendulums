%% RK 4 for Controlled system

function [x,x_sim] = RK4_x(A,B,E,x0,x_k,u_k,K)
   x=zeros(size(x_k));
   x(:,1)=x0;
   h=0.01;
   for i=1:size(x_k,2)-1
        f1 = CS(x(:, i), E(:, :, i), A(:, :, i), B(:, :, i),u_k(i),K(:,:,i),x_k(:,i));
        f2 = CS(x(:, i)+f1/2*h, E(:, :, i), A(:, :, i), B(:, :, i),u_k(i),K(:,:,i),x_k(:,i));
        f3 = CS(x(:, i)+f2/2*h, E(:, :, i), A(:, :, i), B(:, :, i),u_k(i),K(:,:,i),x_k(:,i));
        f4 = CS(x(:, i)+f3*h, E(:, :, i), A(:, :, i), B(:, :, i),u_k(i),K(:,:,i),x_k(:,i));
        x(:, i+1) = x(:, i) + h*(f1/6+((f2+f3)/3)+f4/6);
        x_sim(:,i)=x_k(:,i)+h*x(:,i);
   end
   x_sim(:,length(x_k))=x_k(:,end)+h*x(:,end);
   function dx=CS(x,E,A,B,u_k,K,x_k)
        dx = -(E^-1)*B*K*x_k + (E^-1)*B*u_k + (E^-1)*(A+B*K)*x;
   end
end