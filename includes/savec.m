function savec(filename)

folder=fileparts(filename);
if exist(folder,'dir')==0
    mkdir(folder)    
end
