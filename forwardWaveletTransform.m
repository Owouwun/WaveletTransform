filename = 'C:\Users\user1\source\repos\PicoscopeCpp\PicoscopeCpp\output1.txt';
tau = 'C:\Users\user1\source\repos\PicoscopeCpp\PicoscopeCpp\output2.txt';
s = 'C:\Users\user1\source\repos\PicoscopeCpp\PicoscopeCpp\output3.txt';

delimiter = '\t';

formatSpecTau = '%f%[^\n\r]';
formatSpecS = '%f%[^\n\r]';

fileIDTau = fopen(tau,'r');
fileIDS = fopen(s,'r');

dataArrayTau = textscan(fileIDTau, formatSpecTau, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
dataArrayS = textscan(fileIDS, formatSpecS, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

fclose(fileIDTau);
fclose(fileIDS);

axisTau = [dataArrayTau{1}];
axisS = [dataArrayS{1}];

ncol = size(axisTau);
formatSpec = [ repmat('%f',[1,ncol]), '%[^\n\r]' ];
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'CollectOutput', true, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
fclose(fileID);
Z = [dataArray{1:end-1}];

clearvars filename tau s delimiter formatSpecTau formatSpecS formatSpec fileIDTau fileIDS fileID dataArrayTau dataArrayS dataArray ans;

figure;
surf(axisTau,axisS,Z);
shading interp;
view([0 90]);