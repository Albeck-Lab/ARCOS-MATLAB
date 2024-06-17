classdef validateARCOS
	methods (Static)
		function [meanARIs,data] = validate()
			param1 = {"lifetime",3};
			param2 = {"lifetime",5};
			param3 = {"spreads",200};
			param4 = {"spreads",200,"lifetime",5};
			param5 = {"spreads",300};
			param6 = {"spreads",500};
			
			params = {param1,param2,param3,param4,param5,param6};
			
			meanARIs = cell(length(params),1);
			
			repetitions = 10;     %Repetitions per parameter
			
			data = cell(4,length(params));
			
			
			for param = 1:length(params)
				meanARI = zeros(1,repetitions);
				for iter = 1:repetitions
					 
					[~,x,y,bin,gtlbls] = arcos_utils.gensynthV2(params{param}{:}); 
					lbls = arcos_core(x,y,bin);


					%[spreadTable,x,y,bin,gtlbls] = gensynthV2(params{param}{:});
					%[lbls,~,optionalOut] = arcos_core(x,y,bin);
			
					AR = zeros(1,size(gtlbls,2));
					for i = 1:size(gtlbls,2)
						AR(i) = RandIndexFS(gtlbls(:,i),lbls(:,i));
					end
					meanARI(iter) = mean(AR);
			
			
					%{
					clust_by_time = optionalOut{6};
					clust_by_id = arcos_utils.reformat(clust_by_time);
					raw_data{1}.data.XCoord = x;
					raw_data{1}.data.YCoord = y;
					clust_by_id = arcos_analysis.analyze({clust_by_id},raw_data);
					clust_by_id = clust_by_id{1};
					clear raw_data
					%}
			
				end
				meanARIs{param} = meanARI;
				data{1,param} = gtlbls;
				data{2,param} = lbls;
				data{3,param} = x;
				data{4,param} = y;
			end
		
		end


		function plotSpreads(lbls,x,y,gif)
			clr = hsv(max(lbls,[],'all')+1);
			clr = clr(randperm(length(clr)),:);
			clr(1,:) = 0.95;
		
			set(0,"DefaultFigureWindowStyle","docked")
			for frameToPlot = 1:size(x,2)
    			figure(frameToPlot)
    			gscatter(...
					x(:,frameToPlot),...
					y(:,frameToPlot),...
					lbls(:,frameToPlot),...
					clr(unique(lbls(:,frameToPlot)+1),:),...
					".",20)
    			saveas(gcf,append(pwd,'/',num2str(frameToPlot,'%03.f'),'.png'))
    			close gcf
			end
		
		
			if gif
				files = dir(append(pwd,'/*.png')); %Get list of filenames
				filename = append(pwd,'/','gif.gif'); % Specify the output file name
				for f = 1:length(files) %Loop through files list
    				im = imread(append(files(f).folder,'/',files(f).name)); %Read image into memory
    				[A,map] = rgb2ind(im,256); %Convert from rgb
    				if f == 1
        				imwrite(A,map,filename,'gif','LoopCount',Inf); %Initialize gif
    				else
        				imwrite(A,map,filename,'gif','WriteMode','append'); %Append gif
    				end
    				hold off
				end
				close gcf
			end
		end
	end
end


