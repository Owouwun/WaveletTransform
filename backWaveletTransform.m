filename = 'C:\Users\user1\source\repos\PicoscopeCpp\PicoscopeCpp\output4.txt';
time = 'C:\Users\user1\source\repos\PicoscopeCpp\PicoscopeCpp\output2.txt';
s = 'C:\Users\user1\source\repos\PicoscopeCpp\PicoscopeCpp\output3.txt';

delimiter = '\t';

formatSpec = '%f%[^\n\r]';

fileIDTime = fopen(time,'r');
fileIDS = fopen(s,'r');
fileIDZ = fopen(filename,'r');

dataArrayTime = textscan(fileIDTime, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
dataArrayS = textscan(fileIDS, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
dataArrayZ = textscan(fileIDZ, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

fclose(fileIDTime);
fclose(fileIDS);
fclose(fileIDZ);

axisTime = [dataArrayTime{1}];
axisS = [dataArrayS{1}];
Z = [dataArrayZ{1}];

plot(axisTime,Z);
title('From ' + string(axisS(1)) + ' to ' + axisS(end) + 'Hz');

clearvars filename time delimiter formatSpecAxis ncol s formatSpec fileIDTime fileIDS fileID dataArrayTime dataArrayS dataArray ans;