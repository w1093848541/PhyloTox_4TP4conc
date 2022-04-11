Path = 'C:\Users\10938\Desktop\Jacobson Lab\2022.03\4tp4conc';%Change to the right path before run
File = dir(fullfile(Path,'*tsv'));
FileNames = {File.name};
PeakList = readtable('NAME.txt');
Mass_List = table2array(PeakList(1:end,1));
n = length(Mass_List);
FileNumber = 1;
Filename = FileNames(FileNumber);
Fname = cell2mat(Filename);
fid = fopen(Fname);
GetColumn = readtable(Fname,'FileType','text');
NumC = size(GetColumn,2);
result = zeros(n+1,NumC-7);
DataSet = textscan(fid,repmat('%s',1,NumC));
fclose(fid);
stop = length(DataSet{1,1,1});
t = 1;
for t = 1:n
    result(t+1,1)= Mass_List (t);
end
s = 1;

for s = 1:n
    mass = Mass_List(s);
    i = 7;
    for i = 7:stop
        data_mass1 = DataSet{1,1}{i,1};
        data_mass = str2num(data_mass1);
         diff = mass - data_mass;
        ppm = abs(diff/mass);
        if ppm <= 0.000002
            m = 2;
            for m = 2:NumC-7
                result(s+1,m) = str2num(DataSet{1,m+7}{i,1});
            end
            end
        end
end
result1 = num2cell(result);
s = 2;
for s = 2:NumC-7
    result1{1,s} = DataSet{1,s}{4,1};
end
result2 = result1;
QCs=0;
E33=0;
E32=0;
E31=0;
E23=0;
E22=0;
E21=0;
E13=0;
E12=0;
E11=0;
CT3=0;
CT2=0;
CT1=0;
CT0=0;
s = 2;
%% 
for s = 2:NumC-7
    if result1{1,s} == 'QCs'
        QCs = QCs +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QCs+1} = result1 {m,s};
        end
    end
end
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E33'
        E33 = E33 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QCs+E33+1} = result1 {m,s};
        end
    end
end
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E32'
        E32 = E32 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QCs+E33+E32+1} = result1 {m,s};
        end
    end
end
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E31'
        E31 = E31 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QCs+E33+E32+E31+1} = result1 {m,s};
        end
    end
end
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E23'
        E23 = E23 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QCs+E33+E32+E31+E23+1} = result1 {m,s};
        end
    end
end     
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E22'
        E22 = E22 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QCs+E33+E32+E31+E23+E22+1} = result1 {m,s};
        end
    end
end
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E21'
        E21 = E21 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QCs+E33+E32+E31+E23+E22+E21+1} = result1 {m,s};
        end
    end
end  
  s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E13'
        E13 = E13 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QCs+E33+E32+E31+E23+E22+E21+E13+1} = result1 {m,s};
        end
    end
end    
  s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E12'
        E12 = E12 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QCs+E33+E32+E31+E23+E22+E21+E13+E12+1} = result1 {m,s};
        end
    end
end 
  s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E11'
        E11 = E11 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QCs+E33+E32+E31+E23+E22+E21+E13+E12+E11+1} = result1 {m,s};
        end
    end
end 
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'CT3'
        CT3 = CT3 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QCs+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+1} = result1 {m,s};
        end
    end
end 
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'CT2'
        CT2 = CT2 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QCs+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+CT2+1} = result1 {m,s};
        end
    end
end 
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'CT1'
        CT1 = CT1 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QCs+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+CT2+CT1+1} = result1 {m,s};
        end
    end
end 
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'CT0'
        CT0 = CT0 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QCs+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+CT2+CT1+CT0+1} = result1 {m,s};
        end
    end
end 
%Eliminating metabolites that have more than 50% of missing values in QC
%samples
%% 
s = 2;
NumToDelete = 1;
ToDelete= [];
for s = 2:n
    testQ = 2;
    for testQ = 1:QCs+1
        Test_QC(1,testQ) = result2{s,testQ};
    end
    numberOfZeros = sum(Test_QC == 0);
    if numberOfZeros > QCs/2
       ToDelete (NumToDelete,1) = s;
       NumToDelete = NumToDelete + 1 ;
    end
end
ToDelete = flip(ToDelete);

NumToDelete = 1;
for NumToDelete = 1:length(ToDelete)
    Deleting = ToDelete(NumToDelete,1);
    result2(Deleting,:) = [];
end

Fname(end-3:end) = [];
nWrite = [Fname '_all.xlsx'];
xlswrite(nWrite,result2);