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
W = 15/16;
for k = 1:itermax
    V_new = V_old;
    for i = 1:N
        for j = 1:N
            f = 7*sin(2*pi*x(i))*cos(3*pi*x(i))*sin(2*pi*y(j))*cos(3*pi*y(j));

            [i_00,i_0,i_1,i_2,j_00,j_0,j_1,j_2] = get_coord(i,j,N);

            V_new(i,j)= W/60*(16*V_old(i_0,j)+16*V_old(i_1,j)+16*V_old(i,j_0)+16*V_old(i,j_1)-12*h.^2*f - ...
                V_old(i_2,j)-V_old(i_00,j)-V_old(i,j_2)-V_old(i,j_00))+ (1-W)*V_old(i,j);
        end
    end
    diff = sqrt(sum(sum((V_new - V_old).^2))) / sqrt(sum(sum(V_new.^2)));
    if diff < tol
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

saveas(gcf,'alinea_c.png','png')
%%

% Este script MATLAB permite importar o 
% ficheiro bin?rio criado pelos programas em C
% e produzir um plot. Pode ser muito ?til para
% ver se as condi??es fronteira est?o corretas,
% para ver se n?o h? erros nas fronteiras dos
% subdom?nios, etc.

nx=100;
ny=nx;

fileID = fopen('..\results_2D_c.bin');
array_MPI = fread(fileID, [ny nx],'double');
array_MPI = rot90(array_MPI,1);
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

saveas(gcf,'MPI_alinea_c.png','png')

disp(['erro: ',num2str(sum(sum((V_new-array_MPI).^2))/(100*100))])
disp(['max: ',num2str(max(abs(V_new-array_MPI),[],'all'))])

function [i00,i0,i1,i2,j00,j0,j1,j2] = get_coord(i,j,N)
    j1= j + 1;
    j0 = j - 1;
    j2 = j + 2;
    j00 = j - 2;

    if (j0 <= 0)
        j0 = j0+N;
        j00 = j00+N;
    elseif (j1 > N)
        j1 = j1-N;
        j2 = j2-N;
    elseif (j00 == 0)
        j00 = N;
    elseif (j2 > N)
        j2 = 1;
    end

    i1= i + 1;
    i0 = i - 1;
    i2 = i + 2;
    i00 = i - 2;

    if (i0 <= 0)
        i0 = i0+N;
        i00 = i00+N;
    elseif (i1 > N)
        i1 = i1-N;
        i2 = i2-N;
    elseif (i00 == 0)
        i00 = N;
    elseif (i2 == N+1)
        i2 = 1;
    end
end
