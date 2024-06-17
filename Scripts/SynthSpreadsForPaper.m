[~,x,y,~,~] = arcos_utils.gensynthV2(); %Generate x and y coords

seed = mean(double(evalc('disp("Daniel Oberbauer")')));
rng(seed);
r = -0.5 + (0.5+0.5)*rand(3600,1); % Generate noise matrix

nx = x+r*10; %Inject noise into XCoords;
ny = y+r*10; %Inject noise into YCoords;

idx = randi(3600,2000,1); %Get vector of random indices to remove from dataset to create holes.
nx(idx,:) = [];
ny(idx,:) = [];

d{1}.data.XCoord = nx; %Load the XCoord data into a cell array
d{1}.data.YCoord = ny; %Load the YCoord data into a cell array
bin = arcos_utils.gensynth(d,1,'seed',"Daniel Oberbauer"); %Generate synthetic spreads using the given data


[clust_by_time,clust_by_id,bin,~,~] = arcos(d,1,"EKAR",'bin',bin); %Use ARCOS to detect synthetic spreads
clust_by_id = arcos_analysis.analyze(clust_by_id,d); %Perform ARCOS analysis to get bounds

arcos_plot.plot(clust_by_time,clust_by_id,bin,d,1,1:100) %Plot the spreads.