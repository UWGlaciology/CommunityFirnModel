%%% 12/9/15: I think htis is useless, just a test file

datadir='/Users/maxstev/documents/grad_school/firn/pire/CFM/CommunityFirnModel/code/firnmodel/CompactionCheck';

aadepthcheck=importdata(strcat(datadir,'/depth.csv'));
aaratecheck=importdata(strcat(datadir,'/comprate.csv'));

figure();
clf;
hold on;
plot(aaratecheck(2,2:end),aadepthcheck(2,2:end-1));
set(gca,'ydir','reverse')