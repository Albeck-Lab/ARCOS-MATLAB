a = synth_clust_by_id;
for i = 1:5
	for ii = 1:size(a{i},1)
		for iii = 1:size(a{i}(ii).data,2)
			a{i}(ii).data(iii).class = 2;
		end
	end
end
