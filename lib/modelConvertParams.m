function R = modelConvertParams(R,model);

if numel(fields(R))~=numel(model); 
	error('R and model must be the same size');
end

for mdl = 1:numel(model)
	for k=1:model(mdl).npar
		fstr = str2func(model(mdl).partransform{k});
		R.(model(mdl).name).parConvert(k,:) = fstr(R.(model(mdl).name).E(k,:));
	end
end
