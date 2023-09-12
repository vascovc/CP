clc
clear all
close all

%% 

nx=100;
ny=nx;

% alinea = 'a';
alinea = 'b';
% alinea = 'c';
% alinea = 'd';

% C = false;
% C = true;

f = [alinea '/results_' alinea '.bin'];

% if (C)
%     f = [alinea '/results_' alinea '_C.bin'];
% end

fileID = fopen(f);

array_MPI = fread(fileID, [ny nx],'double');
fclose(fileID);

fexact = [alinea '/vnewMat.mat'];
load(fexact)

L=1;
x=linspace(-L,L,nx);
y=linspace(-L,L,ny);

% figure
% mesh(x,y,array_MPI)
% xlim([-L L])
% ylim([-L L])
% xlabel('\it{x}')
% ylabel('\it{y}')
% title('array\_MPI (no)')
fprintf("erro: %d\n",getMSE(Vnew,array_MPI));

figure
mesh(x,y,array_MPI')
xlim([-L L])
ylim([-L L])
xlabel('\it{x}')
ylabel('\it{y}')
title('array\_MPI')
fprintf("erro trans: %d\n",getMSE(Vnew,array_MPI'));

i = [alinea,'/img', upper(alinea), '.jpg'];
% if (C)
%     i = [alinea,'/img', upper(alinea), '_C.jpg'];
% end
saveas(gcf,i)

% figure
% mesh(x,y,rot90(array_MPI))
% xlim([-L L])
% ylim([-L L])
% xlabel('\it{x}')
% ylabel('\it{y}')
% title('array\_MPI rot90')
fprintf("erro rot90: %d\n",getMSE(Vnew,rot90(array_MPI)));

%%

figure
mesh(x,y,Vnew)
xlim([-L L])
ylim([-L L])
xlabel('\it{x}')
ylabel('\it{y}')
title('array\_Matlab')


%%

function MSE = getMSE(Vnew,array_MPI)

    N2 = 100*100;
    
    
%     matNorm = Vnew./max(Vnew,[],'all');
%     MPINorm = array_MPI./max(array_MPI,[],'all');

    matNorm = Vnew;
    MPINorm = array_MPI;
    
    MSE = 1/N2 * sum((MPINorm-matNorm).^2,'all');
end




