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

h = 2*L/nx;

N = nx;
V_old = zeros(N,N);
% Vold(1,:) = (1 + x)./4;
% Vold(:,1) = (1 + y)./4;
% Vold(end,:) = (3 + x)./4;
% Vold(:,end) = (3 + y)./4;
V_new = V_old;

for k = 1:itermax
    for i = 1:N
        for j = 1:N
            f = 7*sin(2*pi*x(i))*cos(3*pi*x(i))*sin(2*pi*y(j))*cos(3*pi*y(j));  
            [i_0,i_1,j_0,j_1] = get_coord(i,j,N);
            V_new(i,j)= 0.25*(V_new(i_1,j)+V_new(i_0,j)+V_new(i,j_1)+V_new(i,j_0)-h^2*f);
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

saveas(gcf,'alinea_d.png','png')
%%

% Este script MATLAB permite importar o 
% ficheiro bin?rio criado pelos programas em C
% e produzir um plot. Pode ser muito ?til para
% ver se as condi??es fronteira est?o corretas,
% para ver se n?o h? erros nas fronteiras dos
% subdom?nios, etc.

nx=100;
ny=nx;

fileID = fopen('..\results_2D_d.bin');
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

saveas(gcf,'MPI_alinea_d.png','png')

disp(['erro: ',num2str(sum(sum((V_new-array_MPI).^2))/(100*100))])
disp(['max: ',num2str(max(abs(V_new-array_MPI),[],'all'))])

function [i_0,i_1,j_0,j_1] = get_coord(i,j,N)
    i_1= i + 1;
    i_0 = i - 1; 
    j_0 = j - 1;    
    j_1= j + 1;
    if (i_0 <= 0)
        i_0 = i_0+N;
    elseif (i_1 > N)
        i_1 = i_1-N;
    end
    if (j_0 <= 0)
        j_0 = j_0+N;
    elseif (j_1 > N)
        j_1 = j_1-N;
    end
end
