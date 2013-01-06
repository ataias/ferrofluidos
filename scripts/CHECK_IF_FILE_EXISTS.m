function RETURN = CHECK_IF_FILE_EXISTS(strFileName)

%CHECK IF TXT FILE EXISTS
FILE_EXIST = which(strFileName);
if isempty(FILE_EXIST)
  warning('fileNotFound', 'File was not found: %s', strFileName);
  RETURN = 1;
else 
  RETURN = 0;
end;
% -----------------------
