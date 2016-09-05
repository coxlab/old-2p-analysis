function savec(filename)

if filename(end)~='\'
    filename=[filename '\'];
end
folder=fileparts(filename);
if exist(folder,'dir')==0
    mkdir(folder)    
end
