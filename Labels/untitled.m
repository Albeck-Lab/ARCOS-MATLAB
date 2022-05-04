for i = 1:length(table_daniel.daniel)
	table_daniel.class(i) = str2num(cell2mat(table_daniel.daniel(i)));
end
table_daniel.daniel = [];