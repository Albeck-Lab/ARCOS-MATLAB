classdef arcos_ml
	methods (Static)
		function [out,user] = label(clust_by_id,xy)
			disp("Please follow all prompts. If you enter a response unintentionally, please note the cluster # and time #");
			disp("A docked figure will open. Your responses will be in reference to this figure.");
			user = input('What is your name? \n',"s");
			if ischar(user); user = {user}; end % assert user is cell
			prompt = 'Which of these best describes the cluster? 1 = donut, 2 = disk \n';
			fh = figure(); %Figure handle
			set(fh,'WindowStyle','docked') %Set figure window behavior
			for i = 1:numel(xy)
				well = xy(i);
				welldata = clust_by_id{well};
				for clust = 1:size(welldata,1)
					id = welldata(clust).cid;
					clustdata = welldata(clust).data;
					for time = 1:size(clustdata,2)
						clf
						XYCoord = clustdata(time).XYCoord;
						bounds = clustdata(time).bounds;
						plot(XYCoord(:,1),XYCoord(:,2), 'o','MarkerEdgeColor','b','MarkerSize', 3, 'LineStyle', 'none')
						title(append('Cluster: ',int2str(id)," Time: ",int2str(time)));
						xlim([0 1500])
						ylim([0 1500])
						axis square
						hold on
						plot(XYCoord(bounds,1),XYCoord(bounds,2),'r');
						set(gca,'ydir','reverse') %Reverse Y axis (image origin is top left)
						shg
						c = getinput(prompt,[1,2]);
						clust_by_id{well}(clust).data(time).(user{1}) = c;
					end
				end
			end
			close(gcf)
			out = clust_by_id;
		end
		function out = classify(clust_by_id)
% 			table = [];
% 			for i = 1:size(clust_by_id,2)
% 				for ii = 1:size(clust_by_id{i},1)
% 					new = struct2table(clust_by_id{i}(ii).data);
% 					new = removevars(new,{'time','XYCoord','id','bounds','inactive','rocarea','roccount'});
% 					table = vertcat(table, new); %#ok<AGROW> 
% 				end
% 			end
% 			rng("default")
% 			%c = cvpartition(table.id,"Holdout",0.2);
% 			c = cvpartition(table.area,"LeaveOut");
% 			trainingIndices = training(c);
% 			testIndices = test(c);
% 			tableTrain = table(trainingIndices,:);
% 			tableTest = table(testIndices,:);
% 			net = fitcnet(tableTrain,"id");
% 			%testAccuracy = 1 - loss(net,tableTest,"id","LossFun","classiferror");
% 			%confusionchart(tableTest.id,predict(net,tableTest));
% 			%exportONNXNetwork(net,'filename.onnx');
			out = [];
		end
	end
end
function result = getinput(prompt,options)
	a = input(prompt,"s");
	if isempty(a)
		warning('Answer cannot be empty')
		getinput(prompt,options)
	elseif ismember(a,string(options))
		result = a;
	else
		warning(append('Invalid response. Options are ', int2str(options)));
		getinput(prompt,options)
	end
end
