[X,Y]=meshgrid(-18:18,-18:18);
sigma=2;
a=300*exp(-(X.^2+Y.^2)/sigma^2);
b=eye(size(a));
a=conv2(a,b);
a=a+rand(size(a));
[Hxx,Hxy,Hyy,lambda1,lambda2,phi]=hessian(a);

threshold=0.3;

subplot(2,3,1)
imagesc(a); colorbar
title('a')

subplot(2,3,2)
imagesc(phi); colorbar
title('\phi')

subplot(2,3,4)
imagesc(lambda1); colorbar
title('\lambda_1')

subplot(2,3,5)
imagesc(lambda2); colorbar
title('\lambda_2')

subplot(2,3,3)
imagesc(lambda1./(lambda1 - lambda2)); colorbar
title('\lambda_1/(\lambda_1-\lambda_2)')

subplot(2,3,6)
imagesc(lambda2<0 & (lambda1./(lambda1 - lambda2))<threshold); colorbar
title(sprintf('\\lambda_1/(\\lambda_1-\\lambda_2)<%g & \\lambda_2<0',threshold))