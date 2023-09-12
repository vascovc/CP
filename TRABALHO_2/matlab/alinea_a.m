close all
clear all

%% 
itermax = 5e5;
tol = 1e-6;

nx = 100;
ny = nx;
L = 1;
x = linspace(-L,L,nx);
y = linspace(-L,L,ny);

h = 2*L/(nx-1);

N = nx;
V_old = zeros(N,N);

V_old(1,:) = (1 + x)./4;
V_old(:,1) = (1 + y)./4;
V_old(end,:) = (3 + x)./4;
V_old(:,end) = (3 + y)./4;
V_new = V_old;

for k = 1:itermax
    for i = 2:N-1
        for j = 2:N-1
            f = 7*sin(2*pi*x(i))*cos(3*pi*x(i))*sin(2*pi*y(j))*cos(3*pi*y(j));  
            V_new(i,j)= 0.25*(V_old(i+1,j)+V_old(i-1,j)+V_old(i,j+1)+V_old(i,j-1)-h^2*f);
        end
    end
    
    diff = sum(sum((V_new - V_old).^2))/ sum(sum(V_new.^2));
    if sqrt(diff) < tol
        break;
    end
	V_old=V_new;
end 

figure;
mesh(x,y,V_new)
xlim([-L L])
ylim([-L L])
xlabel('x')
ylabel('y')
title('array\_MATLAB')

saveas(gcf,'alinea_a.png','png')
%%

% Este script MATLAB permite importar o 
% ficheiro bin?rio criado pelos programas em C
% e produzir um plot. Pode ser muito ?til para
% ver se as condi??es fronteira est?o corretas,
% para ver se n?o h? erros nas fronteiras dos
% subdom?nios, etc.

nx=100;
ny=nx;

fileID = fopen('..\results_2D_a.bin');
array_MPI = fread(fileID, [ny nx],'double');
array_MPI = rot90(array_MPI,2);
fclose(fileID);

L=1;
x=linspace(-L,L,nx);
y=linspace(L,-L,ny);

% Como para as figuras do MATLAB,
% o avan?o numa linha ? um aumento de x 
% (o MATLAB ? em "column-major"),
% para a figura ficar consistente com o 
% programa em C (o nosso programa C escreve 
% em "row-major"), tem que se transpor a matriz.

figure
mesh(x,y,array_MPI)
xlim([-L L])
ylim([-L L])
xlabel('\it{x}')
ylabel('\it{y}')
title('array\_MPI')

saveas(gcf,'MPI_alinea_a.png','png')

disp(['erro: ',num2str(sum(sum((V_new-array_MPI).^2))/(100*100))])
disp(['max: ',num2str(max(abs(V_new-array_MPI),[],'all'))])