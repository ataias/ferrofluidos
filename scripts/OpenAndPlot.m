% @file OpenAndPlot.m
% @author Ataias Pereira Reis
% OPEN A MATRIX OF DOUBLES IN A BINARY FILE AND PLOTS

function [RETURN, MATRIX] = OpenAndPlot(fileName, nMatrixOrder)
%SELECIONE AS CONDIÇÕES DE ARQUIVO;
FILE_TYPE_BINARY = 1;
FILE_TYPE_ASCII = ~FILE_TYPE_BINARY;
% ----------------

%Check if file exists, if not, the program return 1 and doesn't do anything
%RETURN = CHECK_IF_FILE_EXISTS(fileName);
MATRIX = 0;
RETURN = 0;
if RETURN ==0
% STEP ESPACIAL DX=DY
dDeltaY = 1/(nMatrixOrder-1);
%----------------

% ABRIR ARQUIVO
fileID = fopen(fileName);
Z = fread(fileID, [nMatrixOrder nMatrixOrder], 'double')';
fclose(fileID);
MATRIX = Z;
% --------------

% MALHA XY
dVectorY = 0:dDeltaY:1;
dVectorX = 0:dDeltaY:1;
[X,Y]=meshgrid(dVectorX,dVectorY);
% ----------

%PLOTAR GRÁFICO
figure;
mesh(X,Y,Z);
xlabel('X');
ylabel('Y');
zlabel('Z');
%-----------
end;