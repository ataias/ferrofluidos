% @file ConverterMatrixTextoBinario.m
% @author Ataias Pereira Reis
% @description FUNCTION TO CONVERT A MATRIX IN AN ASCII TXT FILE TO A BINARY FILE
% RETURN = 0 for SUCCESS
% RETURN = 1 for FAILURE

function RETURN = convertMatrixASCIItoBINARY(strFileNameTXT, strFileNameBINARY, nMatrixOrder)

%CHECK IF TXT FILE EXISTS
FILE_EXIST = which(strFileNameTXT);
if isempty(FILE_EXIST)
  warning('fileNotFound', 'File was not found: %s', strFileNameTXT);
  RETURN = 1;
else 
  RETURN=0;
end;
% -----------------------

if RETURN~=1
% READ MATRIX FROM TXT FILE
ReadMatrix = reshape(textread(strFileNameTXT,'%f'),nMatrixOrder,nMatrixOrder);
% -------------------------

% CONVERT AND SAVE TO BINARY FILE
nFileID = fopen(strFileNameBINARY, 'w');
fwrite(nFileID, ReadMatrix', 'double');
fclose(nFileID);
% --------------------------
end;

