Cp=1.0*1000; % J/(kg K)
Tneed=238;
Tair=208:Tneed;
NumHoles=11;

boxsize=0.15^3; % m^3, size of the box enclosing the instrument
rhoair=1.223; %kg/m^3
mair=boxsize*rhoair; %mass of air in box

Q=mair*Cp*(Tneed-Tair); %heat needed to raise temp inside box to Tneed
Qmax=Q(1); %maximum heat needed, i.e. probable coldest to Tneed (-35?)

HeaterSize=10; %Watts
HeatTime=Qmax/HeaterSize; %time needed from heater to heat up to Tneed, seconds;

sysV=12; %Volts
HeatAmp=HeaterSize/sysV; %Amps that the heater draws

HeatAmpHours=HeatTime*HeatAmp/3600; %Amphours for 1 box heating.

HeatAmpHours_yr=HeatAmpHours*365; %Assume 1 measurement per day for the year --> amp hours for the year for heating one box

HeatAmpHours_all=HeatAmpHours_yr*NumHoles; %Total amp hours needed to heat all boxes, 1x per day for the year.

% we are spec'ed for 2W continuous power
% with a 12V system, that is 0.16 Amps continuous
% 8760 hours in a year, so spec'ed for 1401 Amp Hours