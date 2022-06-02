function arcos_id_apoptosis(clust_by_id,xy,ch,nd2path)
	handle=3;
	dao = ImageAccess(nd2path);
	
	for ixy = 1:numel(xy)
		well = xy(ixy);
		spreads = clust_by_id{well};
		for spr = 1:size(spreads,1)
			spread = spreads(spr);
			if spread.t_start-handle>0
				first_frame = spread.t_start-handle;
			else
				first_frame = spread.t_start;
			end
			if spread.t_start+handle>dao.imax.t
				last_frame = dao.imax.t;
			else
				last_frame = spread.t_start+handle;
			end
			frames = first_frame:last_frame;
			for iframe = frames(1):frames(end)
				clf
				CurrFrame = dao.GetFrame([iframe,ch,well,1]);
				CurrFrame = uint16(CurrFrame);
            	adjust = stretchlim(CurrFrame);
            	CurrFrame = imadjust(CurrFrame,adjust,[]);
            	CurrFrame = double(CurrFrame);
            	CurrFrame = uint8(round(CurrFrame/256));
            	CurrFrame = repmat(CurrFrame, [1, 1, 3]);
            	imshow(CurrFrame);
				hold on
				plot(1,1)
				title(append('XY ',string(well),' Spread ',string(spr),' Time: ',string(iframe)))
				spreadData = [];
				for n = 1:size(spread.data,2)
					if iframe == spread.data(n).time
						spreadData = spread.data(n);
						break
					end
				end
				if ~isempty(spreadData)
					bounds = spreadData.bounds;
					points = spreadData.XYCoord;
					hold on
					plot(points(bounds,1),points(bounds,2),'g','LineWidth',1.5)
				end
				im = getframe(gcf);
				
				[A,map] = rgb2ind(im.cdata,256); %Convert from rgb
				filename = append('XY ',string(well),' Spread ',string(spr),'.gif');
				if iframe == frames(1)
                    imwrite(A,map,filename,'gif','LoopCount',Inf); %Initialize gif
                else
                    imwrite(A,map,filename,'gif','WriteMode','append'); %Append gif
				end
				hold off
			end
		end %End spread loop
	end %End XY loop
end