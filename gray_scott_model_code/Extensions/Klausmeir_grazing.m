clear
clc
close all

Nx = 30;
Ny = 30;
%Initial Conditions
w0 = 1 + 0.01*randn(Nx,Ny);
b0 = 1 + 0.01*randn(Nx,Ny);

%Flatten a blank vector to store solutions
u0 = [w0(:);b0(:)];

[t,u] = ode15s(@klausmeir,[0 500],u0);

% --- Plot final state ---
% assume Nx and Ny defined in the main script (30)
N = Nx*Ny;
tf_idx = size(u,1);        % final time index
uf = u(tf_idx,:)';         % column vector

wf = reshape(uf(1:N), Nx, Ny);
bf = reshape(uf(N+1:end), Nx, Ny);

% --- Animate solutions over time ---
figure;

subplot(1,2,1);
im1 = imagesc(zeros(Nx,Ny));
axis image off
colorbar
title('w')

subplot(1,2,2);
im2 = imagesc(zeros(Nx,Ny));
axis image off
colorbar
title('b')

frameIndices = round(linspace(1,size(u,1),100));

for k = frameIndices

    uk = u(k,:)';

    wk = reshape(uk(1:N),Nx,Ny);
    bk = reshape(uk(N+1:end),Nx,Ny);

    set(im1,'CData',wk);
    set(im2,'CData',bk);

    sgtitle(sprintf('t = %.3f',t(k)));

    drawnow
end
function dudt = klausmeir(t,u)

%Parameters
a = 2;
nu = 5;
beta = 0.540;
Dw = 1;
Db = 0.01;
g = 0.4;
k=1;


%Create a Spatial Grid
Nx = 30;
Ny = 30;
L = 10;
dx = L/(Nx-1);
dy = L/(Ny-1);
dt = 0.001;

N = Nx*Ny;
w = reshape(u(1:N),Nx,Ny);
b = reshape(u(N+1:end),Nx,Ny); %Reshape into matrices 

w(1,:)   = w(2,:);
w(end,:) = w(end-1,:);
w(:,1)   = w(:,2);
w(:,end) = w(:,end-1);

b(1,:)   = b(2,:);
b(end,:) = b(end-1,:);
b(:,1)   = b(:,2);
b(:,end) = b(:,end-1);

%Initialise Gradients and Laplacians
gradwx = zeros(Nx,Ny);
gradwy = zeros(Nx,Ny);
lapw = zeros(Nx, Ny);
lapb = zeros(Nx, Ny);

%Initialise solution vectors
dw = zeros(Nx,Ny);
db = zeros(Nx,Ny);

for i = 2:Nx-1
        for j = 2:Ny-1

        %laplacian of w
        lapw(i,j) = (w(i+1,j)-2*w(i,j)+w(i-1,j))/(dx)^2 + (w(i,j+1)-2*w(i,j)+w(i,j-1))/(dy)^2 ;

        %Laplacian of b
        lapb(i,j) = (b(i+1,j)-2*b(i,j)+b(i-1,j))/(dx)^2  + (b(i,j+1)-2*b(i,j)+b(i,j-1))/(dy)^2 ;

        %nabla of w
        gradwx(i,j) = (w(i+1,j) - w(i-1,j))/(2*dx);
        gradwy(i,j) = (w(i,j+1) - w(i,j-1))/(2*dy);

        %Compute RHS of PDES
        dw(i,j) = a - w(i,j) - w(i,j)*b(i,j)^2 + nu*gradwx(i,j) + Dw*lapw(i,j);
        db(i,j) = w(i,j)*b(i,j)^2 - beta*b(i,j) + Db*lapb(i,j) - (g*b(i,j))/(k+b(i,j));
        end
end

dw(1,:)   = dw(2,:);
dw(end,:) = dw(end-1,:);
dw(:,1)   = dw(:,2);
dw(:,end) = dw(:,end-1);

db(1,:)   = db(2,:);
db(end,:) = db(end-1,:);
db(:,1)   = db(:,2);
db(:,end) = db(:,end-1);

dudt = [dw(:); db(:)];
end
