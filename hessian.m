function [Hxx,Hxy,Hyy,lambda1,lambda2]=hessian(a)
% compute the hessian matrix of an equispaced image matrix, and its 
%  eigenvalues, by central finite differences
% Neglect edges for simplicity - h, lambda1,2 are two rows and two columns
%  smaller than a
% Omit for simplicity the scaling 1/4dxdy

[ny,nx]=size(a);

Hxx = a(1:ny-2,3:nx) - 2*a(1:ny-2,2:nx-1) + a(1:ny-2,1:nx-2);
Hyy = a(3:ny,1:nx-2) - 2*a(2:ny-1,1:nx-2) + a(1:ny-2,1:nx-2);
Hxy = a(3:ny,3:nx) - a(3:ny,1:nx-2) - a(1:ny-2,3:nx) + a(1:ny-2,1:nx-2);

D = sqrt((Hxx-Hyy).^2 + 4*Hxy.^2);
lambda1 = (Hxx + Hyy + D)/2;
lambda2 = (Hxx + Hyy - D)/2;
