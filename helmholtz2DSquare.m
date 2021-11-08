% 2D Helmholtz
clear; clc % Clear workspace and command window

%Inputs

Nx=100; % 50 grid points in x and y 
Ny=100;
Lx=5; % Domain size
Ly=5;

% Grid - Eg. 
% Range from 0 to Lx with Nx grid points

x = linspace(0,Lx,Nx); % Mesh
y = linspace(0,Ly,Ny);

dx = x(2) - x(1); % Mesh Size
dy = y(2) - y(1);

% Initialise matrices

N = Nx*Ny; % no. of unknowns

M = zeros(N,N); % N rows, N columns

% Interior points
for i=2:Nx-1 % Loop over x direction skipping 1st and last grid points
    for j=2:Ny-1 % Loop over y direction skipping first and last grid points
        n=i+(j-1)*Nx; % Convert ij grid point to nth grid point (n refers to the nth entry of the column vector in which the matrix acts upon)
        M(n,n)=-2*((1/(dx)^2) + (1/(dy)^2)); % Main diagonal
        M(n,n-1)=1/((dx)^2); % Off diagonal to the left ie. point i-1, j
        M(n,n+1)=1/((dx)^2); % Off diagonal to the right ie. point i+1, j
        M(n,n-Nx)=1/((dy)^2); % Far off diagonal to the left ie. point i, j-1
        M(n,n+Nx)=1/((dy)^2); % Far off diagonal to the left ie. point i, j+1
    end
end

% Boundary Conditions

% Left BC: 
i=1;
for j=1:Ny
    n=i+(j-1)*Nx; % nth row for this ij
    M(n,n)=1; % Main diagonal
    M(n,n+1)=-1; %Off diagonal to the right ie. point i+1, j
end
     
% Right BC: 
i=Nx-1;
for j=1:Ny
    n=i+(j-1)*Nx; % nth row for this ij
    M(n,n)=1; % Main diagonal
    M(n,n-1)=-1;
end

% Bottom BC: 
j=1;
for i=1:Nx
    n=i+(j-1)*Nx; % nth row for this ij
    M(n,n)=1; % Main diagonal
    M(n,n+Nx)=-1;
end

% Top BC: 
j=Ny-1;
for i=1:Nx
    n=i+(j-1)*Nx; % nth row for this ij
    M(n,n)=1; % Main diagonal
    M(n,n-Nx)=-1;    
end

% Solve for z
[z_vec,alpha] = eigs(M);

z = zeros(Nx,Ny);

% Convert vector z back into matrix
for i=1:Nx
     for j=1:Ny
         n=i+(j-1)*Nx;
         z(i,j) = z_vec(n,2);
     end
end

% Plot result

surf(x,y,z') % ' takes the transpose (Matlab thinks about vectors weirdly)
xlabel('x')
ylabel('y')
set(gca, 'Fontsize',16)
