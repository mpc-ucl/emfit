function cleanpath(modelClass)

warning('off','MATLAB:rmpath:DirNotFound');
for k=1:length(modelClass)
	rmpath(genpath(modelClass{k}));
end
