folderPath = 'F:\Sentinel-3-Level1\Lakes\047 Mohavi\L1BS';
L1BS=Load_L1BS_Data_products(folderPath);
MOHAVE='S3AMohave_p3101t95o110sr3000wo80or99.mat';
load(MOHAVE);
dataset=S3AMohave_p3101t95o110sr3000wo80or99;
load('09422500.mat')
ref=USGS_Data{2,2};
USGS_POS=[ 	35.1963472,	-114.5702028];

cycles_L1BS=[];
Raw_elevation_total=dataset.ObjVS.Raw.Sat.Alt.Hi.Ku.Signal-dataset.ObjVS.Raw.Mes.Rng.Tracker.SAR.Ku.Signal;

ck=dataset.ObjVS.Raw.Sat.Alt.Hi.Ku.Cycle';
cycles_alt_bundle=unique(ck);

for i=1:numel(L1BS)
    cycleNumbers = char(regexp(string(L1BS(i)), '_\d{3}_', 'match'));
    tmp=str2num(cycleNumbers(2:4));
    cycles_L1BS=[cycles_L1BS;tmp];
end
