% Code to solve 1D Helmholtz Equation
clear; clc % Clear workspace and command window

% Inputs
Nx=500; % 500 grid points in x
Lx=2*pi; % Domain size

% Grid 
% Range from 0 to Lx with Nx grid points

x = linspace(0,Lx,Nx); % Mesh
h = x(2) - x(1); % Mesh Size

s = 1; % Constant

% Initialise matrices (eqn: Mz = B)

N = Nx; % No. of unknowns

M = zeros(N,N); % N rows, N columns
B = zeros(N,1); % N rows, 1 column

% Left BC: at first point
i=1;
n=i; % nth row for this ij
M(n,n)=1; % Main diagonal
% M(n,n-1)=-1; % Off diagonal to the left 
% M(n,n+1)=-1; % off diagonal to the right
B(n,1)=sin((x(i)*pi)/2); % BC 
     
% Right BC: at Lx
i=Nx;
n=i; % nth row for this ij
M(n,n)=1; % Main diagonal
% M(n,n-1)=-1; % Off diagonal to the left 
% M(n,n+1)=1; % off diagonal to the right
B(n,1)=sin((x(i)*pi)/2); % BC 

% Interior points
for i=2:Nx-1 % Loop over x direction skipping 1st and last grid points
        n=i; % Convert ij grid point to nth grid point (n refers to the nth entry of the column vector in which the matrix acts upon)
        M(n,n)=2-(s^2*h^2); % Main diagonal
        M(n,n-1)=-1; % Off diagonal to the left 
        M(n,n+1)=-1; % off diagonal to the right
        % M(n,n-Nx)=0; % Far off diagonal to the left
        % M(n,n+Nx)=0; % Far off diagonal to the left
        B(n,1) = 0; % Source term 
end


% Solve for z
z_vec=M\B;

% Convert vector z back into a 2D array
z = zeros(1,n);
for i=1:Nx
    n=i;
    z(i) = z_vec(n);
end

% Plot result 
% surf(x,y,z') % ' takes the transpose (Matlab thinks about vectors weirdly)
% xlabel('x')
% set(gca, 'Fontsize', 16)
plot(x,z)


