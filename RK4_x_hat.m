%% RK 4 for Estimated system

function [x_hat,x_hat_sim] = RK4_x_hat(A,B,L,C,E,x_hat_0,x_k,u_k,K,x)   
   x_hat=zeros(size(x_k));
   x_hat(:,1)=x_hat_0;
   h=0.01;
   for i=1:size(x_k,2)-1
        f1 = ES(x_hat(:, i), E(:, :, i), A(:, :, i), B(:, :, i),C(:, :, i),u_k(i),K(:,:,i),L(:,:,i),x_k(:,i),x(:,i));
        f2 = ES(x_hat(:, i)+f1/2*h, E(:, :, i), A(:, :, i), B(:, :, i),C(:, :, i),u_k(i),K(:,:,i),L(:,:,i),x_k(:,i),x(:,i));
        f3 = ES(x_hat(:, i)+f2/2*h, E(:, :, i), A(:, :, i), B(:, :, i),C(:, :, i),u_k(i),K(:,:,i),L(:,:,i),x_k(:,i),x(:,i));
        f4 = ES(x_hat(:, i)+f3*h, E(:, :, i), A(:, :, i), B(:, :, i),C(:, :, i),u_k(i),K(:,:,i),L(:,:,i),x_k(:,i),x(:,i));
        x_hat(:, i+1) = x_hat(:, i) + h*(f1/6+((f2+f3)/3)+f4/6);
        x_hat_sim(:,i)=x_k(:,i)+h*x_hat(:,i);
   end
   x_hat_sim(:,length(x_k))=x_k(:,end)+h*x_hat(:,end);

   function dx_hat=ES(x_hat,E,A,B,C,u_k,K,L,x_k,x)
        dx_hat = (E^-1*A + L*C)*x_hat + (2*L*C - E^-1*B*K)*x_k + (L*C + E^-1*B*K)*x + E^-1*B*u_k;
   end
end