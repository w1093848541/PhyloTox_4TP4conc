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
result3 = readmatrix('PQN_predictedresult.csv');
QC=0;
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

for s = 2:NumC-7
    if result1{1,s} == 'QCs'
        QC = QC +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+1} = result1 {m,s};
        end
    end
end
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E33'
        E33 = E33 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E33+1} = result1 {m,s};
        end
    end
end
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E32'
        E32 = E32 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E33+E32+1} = result1 {m,s};
        end
    end
end
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E31'
        E31 = E31 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E33+E32+E31+1} = result1 {m,s};
        end
    end
end
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E23'
        E23 = E23 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E33+E32+E31+E23+1} = result1 {m,s};
        end
    end
end     
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E22'
        E22 = E22 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E33+E32+E31+E23+E22+1} = result1 {m,s};
        end
    end
end
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E21'
        E21 = E21 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E33+E32+E31+E23+E22+E21+1} = result1 {m,s};
        end
    end
end  
  s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E13'
        E13 = E13 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E33+E32+E31+E23+E22+E21+E13+1} = result1 {m,s};
        end
    end
end    
  s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E12'
        E12 = E12 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E33+E32+E31+E23+E22+E21+E13+E12+1} = result1 {m,s};
        end
    end
end 
  s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'E11'
        E11 = E11 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+1} = result1 {m,s};
        end
    end
end 
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'CT3'
        CT3 = CT3 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+1} = result1 {m,s};
        end
    end
end 
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'CT2'
        CT2 = CT2 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+CT2+1} = result1 {m,s};
        end
    end
end 
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'CT1'
        CT1 = CT1 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+CT2+CT1+1} = result1 {m,s};
        end
    end
end 
s = 2;
for s = 2:NumC-7
    if result1{1,s} == 'CT0'
        CT0 = CT0 +1;
        m = 1;
        for m = 1:n+1
            result2 {m,QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+CT2+CT1+CT0+1} = result1 {m,s};
        end
    end
end 
% Fname(end-3:end) = [];
% nWrite = [Fname '_all.xlsx'];
% xlswrite(nWrite,result2);

%% 


Col = size (result3,2);
Row = size (result3,1);
result3 = num2cell(result3);
ColResult2 = 2;
for ColResult2 = 2:Col
    RowResult2 = 2;
    for RowResult2 = 2: Row+1
        result2 {RowResult2,ColResult2} = result3{RowResult2-1,ColResult2};
    end
end

% Below is the t-test 


resultfort = cell2mat(result3.');
t = 1;



n = length(result3);
s = size(result3,2);
ttestallT3_CV3 = 1;
ttestallT3_CV2 = 1;
ttestallT3_CV1 = 1;
ttestallT3_2V3 = 1;
ttestallT3_2V1 = 1;
ttestallT3_1V3 = 1;
ttestallT2_CV3 = 1;
ttestallT2_CV2 = 1;
ttestallT2_CV1 = 1;
ttestallT2_2V3 = 1;
ttestallT2_2V1 = 1;
ttestallT2_1V3 = 1;
ttestallT1_CV3 = 1;
ttestallT1_CV2 = 1;
ttestallT1_CV1 = 1;
ttestallT1_2V3 = 1;
ttestallT1_2V1 = 1;
ttestallT1_1V3 = 1;

for t = 1:n
    QCt = zeros(QC,1);
    E33t= zeros(E33,1);
    E32t= zeros(E32,1);
    E31t= zeros(E31,1);
    E23t= zeros(E23,1);
    E22t= zeros(E22,1);
    E21t= zeros(E21,1);
    E13t= zeros(E13,1);
    E12t= zeros(E12,1);
    E11t= zeros(E11,1);
    CT3t= zeros(CT3,1);
    CT2t= zeros(CT2,1);
    CT1t= zeros(CT1,1);
    CT0t= zeros(CT0,1);
    
 
    filt = 1;
    for filt = 1:QC+1
        QCt(filt,1) = resultfort(filt,t);
    end
    filt = QC+1;
    for filt = QC+1:QC+E33+1
        E33t(filt,1) = resultfort(filt,t);
    end
    filt = QC+E33+1;
    for filt = QC+E33+1:QC+E33+E32+1
        E32t(filt,1) = resultfort(filt,t);
    end
    filt = QC+E33+E32+1;
    for filt = QC+E33+E32+1:QC+E33+E32+E31+1
        E31t(filt,1) = resultfort(filt,t);
    end
    filt = QC+E33+E32+E31+1;
    for filt = QC+E33+E32+E31+1:QC+E33+E32+E31+E23+1
        E23t(filt,1) = resultfort(filt,t);
    end
    filt = QC+E33+E32+E31+E23+1;
    for filt = QC+E33+E32+E31+E23+1:QC+E33+E32+E31+E23+E22+1
        E22t(filt,1) = resultfort(filt,t);
    end
    filt = QC+E33+E32+E31+E23+E22+1;
    for filt = QC+E33+E32+E31+E23+E22+1:QC+E33+E32+E31+E23+E22+E21+1
        E21t(filt,1) = resultfort(filt,t);
    end
    filt = QC+E33+E32+E31+E23+E22+E21+1;
    for filt = QC+E33+E32+E31+E23+E22+E21+1:QC+E33+E32+E31+E23+E22+E21+E13+1
        E13t(filt,1) = resultfort(filt,t);
    end 
    filt = QC+E33+E32+E31+E23+E22+E21+E13+12;
    for filt = QC+E33+E32+E31+E23+E22+E21+E13+1:QC+E33+E32+E31+E23+E22+E21+E13+E12+1
        E12t(filt,1) = resultfort(filt,t);
    end
    filt = QC+E33+E32+E31+E23+E22+E21+E13+E12+1;
    for filt = QC+E33+E32+E31+E23+E22+E21+E13+E12+1:QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+1
        E11t(filt,1) = resultfort(filt,t);
    end
    filt = QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+1;
    for filt = QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+1:QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+1
        CT3t(filt,1) = resultfort(filt,t);
    end
    filt = QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+1;
    for filt = QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+1: QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+CT2+1
        CT2t(filt,1) = resultfort(filt,t);
    end
    filt = QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+CT2+1;
    for filt = QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+CT2+1: QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+CT2+CT1+1
        CT1t(filt,1) = resultfort(filt,t);
    end
    filt = QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+CT2+CT1+1;
    for filt = QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+CT2+CT1+1: QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+CT2+CT1+CT0+1
        CT0t(filt,1) = resultfort(filt,t);
    end
    
    
    QCt(QCt==0) = [];
    E33t(E33t==0) = [];
    E32t(E32t==0) = [];
    E31t(E31t==0) = [];
    E23t(E23t==0) = [];
    E22t(E22t==0) = [];
    E21t(E21t==0) = [];
    E13t(E13t==0) = [];
    E12t(E12t==0) = [];
    E11t(E11t==0) = [];
    CT3t(CT3t==0) = [];
    CT2t(CT2t==0) = [];
    CT1t(CT1t==0) = [];
    CT0t(CT0t==0) = [];
  
    
    [h,p] = ttest2(E33t,E23t,'Vartype','unequal');
    meanE33t = mean(E33t(:));
    meanE23t = mean(E23t(:));
    log2FC = log2(meanE33t/meanE23t);
       ts = 1;
       for ts = 1:s
       ttest2resultT3_2V3(ttestallT3_2V3,ts) = cell2mat(result3(t,ts));
       ttest2resultT3_2V3(ttestallT3_2V3,s+1) = p;
       ttest2resultT3_2V3(ttestallT3_2V3,s+2) = log2FC;
       end
       ttestallT3_2V3 = ttestallT3_2V3 +1;
     
       
     [h,p] = ttest2(E33t,E13t,'Vartype','unequal');
     meanE33t = mean(E33t(:));
     meanE13t = mean(E13t(:));
     log2FC = log2(meanE33t/meanE13t);
       ts = 1;
       for ts = 1:s
       ttest2resultT3_1V3(ttestallT3_1V3,ts) = cell2mat(result3(t,ts));
       ttest2resultT3_1V3(ttestallT3_1V3,s+1) = p;
       ttest2resultT3_1V3(ttestallT3_1V3,s+2) = log2FC;
       end
       ttestallT3_1V3 = ttestallT3_1V3 +1;
       
        [h,p] = ttest2(E33t,CT3t,'Vartype','unequal');
    meanE33t = mean(E33t(:));
    meanCT3t = mean(CT3t(:));
    log2FC = log2(meanE33t/meanCT3t);
       ts = 1;
       for ts = 1:s
       ttest2resultT3_CV3(ttestallT3_CV3,ts) = cell2mat(result3(t,ts));
       ttest2resultT3_CV3(ttestallT3_CV3,s+1) = p;
       ttest2resultT3_CV3(ttestallT3_CV3,s+2) = log2FC;
       end
       ttestallT3_CV3 = ttestallT3_CV3 +1;
       
        [h,p] = ttest2(E23t,E13t,'Vartype','unequal');
    meanE23t = mean(E23t(:));
    meanE13t = mean(E13t(:));
    log2FC = log2(meanE23t/meanE13t);
       ts = 1;
       for ts = 1:s
       ttest2resultT3_1V2(ttestallT3_2V1,ts) = cell2mat(result3(t,ts));
       ttest2resultT3_1V2(ttestallT3_2V1,s+1) = p;
       ttest2resultT3_1V2(ttestallT3_2V1,s+2) = log2FC;
       end
       ttestallT3_2V1 = ttestallT3_2V1 +1;
       
       
        [h,p] = ttest2(E23t,CT3t,'Vartype','unequal');
    meanE23t = mean(E23t(:));
    meanCT3t = mean(CT3t(:));
    log2FC = log2(meanE23t/meanCT3t);
       ts = 1;
       for ts = 1:s
       ttest2resultT3_CV2(ttestallT3_CV2,ts) = cell2mat(result3(t,ts));
       ttest2resultT3_CV2(ttestallT3_CV2,s+1) = p;
       ttest2resultT3_CV2(ttestallT3_CV2,s+2) = log2FC;
       end
       ttestallT3_CV2 = ttestallT3_CV2 +1;
       
       
        [h,p] = ttest2(E13t,CT3t,'Vartype','unequal');
    meanE13t = mean(E13t(:));
    meanCT3t = mean(CT3t(:));
    log2FC = log2(meanE13t/meanCT3t);
       ts = 1;
       for ts = 1:s
       ttest2resultT3_CV1(ttestallT3_CV1,ts) = cell2mat(result3(t,ts));
       ttest2resultT3_CV1(ttestallT3_CV1,s+1) = p;
       ttest2resultT3_CV1(ttestallT3_CV1,s+2) = log2FC;
       end
       ttestallT3_CV1 = ttestallT3_CV1 +1;
       
       
       
       
       
       
    [h,p] = ttest2(E32t,E22t,'Vartype','unequal');
    meanE32t = mean(E32t(:));
    meanE22t = mean(E22t(:));
    log2FC = log2(meanE32t/meanE22t);
       ts = 1;
       for ts = 1:s
       ttest2resultT2_2V3(ttestallT2_2V3,ts) = cell2mat(result3(t,ts));
       ttest2resultT2_2V3(ttestallT2_2V3,s+1) = p;
       ttest2resultT2_2V3(ttestallT2_2V3,s+2) = log2FC;
       end
       ttestallT2_2V3 = ttestallT2_2V3 +1;
     
       
     [h,p] = ttest2(E32t,E12t,'Vartype','unequal');
     meanE32t = mean(E32t(:));
     meanE12t = mean(E12t(:));
     log2FC = log2(meanE32t/meanE12t);
       ts = 1;
       for ts = 1:s
       ttest2resultT2_1V3(ttestallT2_1V3,ts) = cell2mat(result3(t,ts));
       ttest2resultT2_1V3(ttestallT2_1V3,s+1) = p;
       ttest2resultT2_1V3(ttestallT2_1V3,s+2) = log2FC;
       end
       ttestallT2_1V3 = ttestallT2_1V3 +1;
       
        [h,p] = ttest2(E32t,CT2t,'Vartype','unequal');
    meanE32t = mean(E32t(:));
    meanCT2t = mean(CT2t(:));
    log2FC = log2(meanE32t/meanCT2t);
       ts = 1;
       for ts = 1:s
       ttest2resultT2_CV3(ttestallT2_CV3,ts) = cell2mat(result3(t,ts));
       ttest2resultT2_CV3(ttestallT2_CV3,s+1) = p;
       ttest2resultT2_CV3(ttestallT2_CV3,s+2) = log2FC;
       end
       ttestallT2_CV3 = ttestallT2_CV3 +1;
       
        [h,p] = ttest2(E22t,E12t,'Vartype','unequal');
    meanE22t = mean(E22t(:));
    meanE12t = mean(E12t(:));
    log2FC = log2(meanE22t/meanE12t);
       ts = 1;
       for ts = 1:s
       ttest2resultT2_1V2(ttestallT2_2V1,ts) = cell2mat(result3(t,ts));
       ttest2resultT2_1V2(ttestallT2_2V1,s+1) = p;
       ttest2resultT2_1V2(ttestallT2_2V1,s+2) = log2FC;
       end
       ttestallT2_2V1 = ttestallT2_2V1 +1;
       
       
        [h,p] = ttest2(E22t,CT2t,'Vartype','unequal');
    meanE22t = mean(E22t(:));
    meanCT2t = mean(CT2t(:));
    log2FC = log2(meanE22t/meanCT2t);
       ts = 1;
       for ts = 1:s
       ttest2resultT2_CV2(ttestallT2_CV2,ts) = cell2mat(result3(t,ts));
       ttest2resultT2_CV2(ttestallT2_CV2,s+1) = p;
       ttest2resultT2_CV2(ttestallT2_CV2,s+2) = log2FC;
       end
       ttestallT2_CV2 = ttestallT2_CV2 +1;
       
       
        [h,p] = ttest2(E12t,CT2t,'Vartype','unequal');
    meanE12t = mean(E12t(:));
    meanCT2t = mean(CT2t(:));
    log2FC = log2(meanE12t/meanCT2t);
       ts = 1;
       for ts = 1:s
       ttest2resultT2_CV1(ttestallT2_CV1,ts) = cell2mat(result3(t,ts));
       ttest2resultT2_CV1(ttestallT2_CV1,s+1) = p;
       ttest2resultT2_CV1(ttestallT2_CV1,s+2) = log2FC;
       end
       ttestallT2_CV1 = ttestallT2_CV1 +1;
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
    [h,p] = ttest2(E31t,E21t,'Vartype','unequal');
    meanE31t = mean(E31t(:));
    meanE21t = mean(E21t(:));
    log2FC = log2(meanE31t/meanE21t);
       ts = 1;
       for ts = 1:s
       ttest2resultT1_2V3(ttestallT1_2V3,ts) = cell2mat(result3(t,ts));
       ttest2resultT1_2V3(ttestallT1_2V3,s+1) = p;
       ttest2resultT1_2V3(ttestallT1_2V3,s+2) = log2FC;
       end
       ttestallT1_2V3 = ttestallT1_2V3 +1;
     
       
     [h,p] = ttest2(E31t,E11t,'Vartype','unequal');
     meanE31t = mean(E31t(:));
     meanE11t = mean(E11t(:));
     log2FC = log2(meanE31t/meanE11t);
       ts = 1;
       for ts = 1:s
       ttest2resultT1_1V3(ttestallT1_1V3,ts) = cell2mat(result3(t,ts));
       ttest2resultT1_1V3(ttestallT1_1V3,s+1) = p;
       ttest2resultT1_1V3(ttestallT1_1V3,s+2) = log2FC;
       end
       ttestallT1_1V3 = ttestallT1_1V3 +1;
       
        [h,p] = ttest2(E31t,CT1t,'Vartype','unequal');
    meanE31t = mean(E31t(:));
    meanCT1t = mean(CT1t(:));
    log2FC = log2(meanE31t/meanCT1t);
       ts = 1;
       for ts = 1:s
       ttest2resultT1_CV3(ttestallT1_CV3,ts) = cell2mat(result3(t,ts));
       ttest2resultT1_CV3(ttestallT1_CV3,s+1) = p;
       ttest2resultT1_CV3(ttestallT1_CV3,s+2) = log2FC;
       end
       ttestallT1_CV3 = ttestallT1_CV3 +1;
       
        [h,p] = ttest2(E21t,E11t,'Vartype','unequal');
    meanE21t = mean(E21t(:));
    meanE11t = mean(E11t(:));
    log2FC = log2(meanE21t/meanE11t);
       ts = 1;
       for ts = 1:s
       ttest2resultT1_1V2(ttestallT1_2V1,ts) = cell2mat(result3(t,ts));
       ttest2resultT1_1V2(ttestallT1_2V1,s+1) = p;
       ttest2resultT1_1V2(ttestallT1_2V1,s+2) = log2FC;
       end
       ttestallT1_2V1 = ttestallT1_2V1 +1;
       
       
        [h,p] = ttest2(E21t,CT1t,'Vartype','unequal');
    meanE21t = mean(E21t(:));
    meanCT1t = mean(CT1t(:));
    log2FC = log2(meanE21t/meanCT1t);
       ts = 1;
       for ts = 1:s
       ttest2resultT1_CV2(ttestallT1_CV2,ts) = cell2mat(result3(t,ts));
       ttest2resultT1_CV2(ttestallT1_CV2,s+1) = p;
       ttest2resultT1_CV2(ttestallT1_CV2,s+2) = log2FC;
       end
       ttestallT1_CV2 = ttestallT1_CV2 +1;
       
       
        [h,p] = ttest2(E11t,CT1t,'Vartype','unequal');
    meanE11t = mean(E11t(:));
    meanCT1t = mean(CT1t(:));
    log2FC = log2(meanE11t/meanCT1t);
       ts = 1;
       for ts = 1:s
       ttest2resultT1_CV1(ttestallT1_CV1,ts) = cell2mat(result3(t,ts));
       ttest2resultT1_CV1(ttestallT1_CV1,s+1) = p;
       ttest2resultT1_CV1(ttestallT1_CV1,s+2) = log2FC;
       end
       ttestallT1_CV1 = ttestallT1_CV1 +1;
       
        
      
 
    
    
       
    
       
    
    
 
       
    
end    

%Output the t-test p-values and log2FC%

nWrite = ['ttest2resultT1_CV1_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT1_CV1);

nWrite = ['ttest2resultT1_CV2_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT1_CV2);

nWrite = ['ttest2resultT1_CV3_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT1_CV3);

nWrite = ['ttest2resultT1_1V2_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT1_1V2);

nWrite = ['ttest2resultT1_1V3_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT1_1V3);

nWrite = ['ttest2resultT1_2V3_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT1_2V3);





nWrite = ['ttest2resultT2_CV1_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT2_CV1);

nWrite = ['ttest2resultT2_CV2_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT2_CV2);

nWrite = ['ttest2resultT2_CV3_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT2_CV3);

nWrite = ['ttest2resultT2_1V2_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT2_1V2);

nWrite = ['ttest2resultT2_1V3_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT2_1V3);

nWrite = ['ttest2resultT2_2V3_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT2_2V3);







nWrite = ['ttest2resultT3_CV1_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT3_CV1);

nWrite = ['ttest2resultT3_CV2_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT3_CV2);

nWrite = ['ttest2resultT3_CV3_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT3_CV3);

nWrite = ['ttest2resultT3_1V2_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT3_1V2);

nWrite = ['ttest2resultT3_1V3_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT3_1V3);

nWrite = ['ttest2resultT3_2V3_ALL.xlsx'];

xlswrite(nWrite,ttest2resultT3_2V3);





%% 

%Opearte the FDR test and output the results with q-values, saved for
%volcano plots


File = dir(fullfile(Path,'*xlsx'));
FileNames = {File.name};


FileCount = size(File,1);
for fc = 1: FileCount
    FDRpass_all = [];
    Filename = FileNames(fc);
    Fname = cell2mat(Filename);
    
    if Fname(1:6) == "ttest2"
    Pvaluesall  =  xlsread(Fname,'Sheet1');
    Pvalues = Pvaluesall(:,Col+1);
    log2FC = Pvaluesall(:,Col+2);
    ids = Pvaluesall(:,1);
    fdr = mafdr(Pvalues);
    [fdr,q,priori,R2] = mafdr(Pvalues,'Method','polynomial');

    qvalue_result4 = [ids,fdr,q,log2FC];

    priori_R24 = [priori,R2];
    Fname(end-8:end) = [];
    Fname(1:12) = [];
    
    for nameid = 1:size(qvalue_result4,1)
        nameid_c = num2str(nameid);
        id_list{nameid} = append('Metabolite_',nameid_c);
    end
    qvalue_result4 = num2cell(qvalue_result4);
    id_list = rot90(id_list);
    id_list = rot90(id_list);
    id_list = rot90(id_list);
    qvalue_result4 = [qvalue_result4,id_list];
    
    id_list = {};
    
    nWrite = [Fname '_qvalues.xlsx'];

    xlswrite(nWrite,qvalue_result4);

    nWrite = [Fname '_priori_R2.xlsx'];

    xlswrite(nWrite,priori_R24);
    
    FDRpass = 1;

    for npass = 1:n
        if cell2mat(qvalue_result4(npass,3)) <= 0.05 
        
        FDRpass_all(FDRpass,:) = cell2mat(result3(npass,:));
      
        FDRpass = FDRpass + 1;
        end
    end
%     p_position = size(FDRpass_all,2);
%     p_position = p_position -1;
%     FDRpass_all(:,p_position) = [];
%     FDRpass_all(:,p_position) = [];
    if size(FDRpass_all) >=1
    nWrite = [Fname '_qpass.xlsx'];

    xlswrite(nWrite,FDRpass_all);
    end
    
    end
    
    qvalue_result4=[];
end






%Collect the metabolites that passed the FDR test

% FDRpass = 1;
% 
% for npass = 1:n
%     if qvalue_result1(npass,3) <= 0.05 || qvalue_result2(npass,3) <= 0.05 ||qvalue_result3(npass,3) <= 0.05 ||qvalue_result4(npass,3) <= 0.05
%        FDRpass_all(FDRpass,:) = ttest2resultT1(FDRpass,:);
%       
%        FDRpass = FDRpass + 1;
%     end
% end
% p_position = size(FDRpass_all,2);
% p_position = p_position -1;
% FDRpass_all(:,p_position) = [];
% FDRpass_all(:,p_position) = [];


% Below is the output for each time point


% FDRpass1 = 1;
% 
% for npass = 1:n
%     if qvalue_result1(npass,3) <= 0.05 
%        FDRpass_1(FDRpass1,:) = ttest2resultT1(FDRpass1,:);
%       
%        FDRpass1 = FDRpass1 + 1;
%     end
% end
% 
% FDRpass2 = 1;
% 
% for npass = 1:n
%     if qvalue_result2(npass,3) <= 0.05
%        FDRpass_2(FDRpass2,:) = ttest2resultT2(FDRpass2,:);
%       
%        FDRpass2 = FDRpass2 + 1;
%     end
% end
% 
% 
% FDRpass3 = 1;
% 
% for npass = 1:n
%     if qvalue_result3(npass,3) <= 0.05
%        FDRpass_3(FDRpass3,:) = ttest2resultT3(FDRpass3,:);
%       
%        FDRpass3 = FDRpass3 + 1;
%     end
% end
% 
% FDRpass4 = 1;
% 
% for npass = 1:n
%     if qvalue_result4(npass,3) <= 0.05
%        FDRpass_4(FDRpass4,:) = ttest2resultT4(FDRpass4,:);
%       
%        FDRpass4 = FDRpass4 + 1;
%     end
% end
% FDRpass_1(:,p_position) = [];
% FDRpass_1(:,p_position) = [];
% 
% FDRpass_2(:,p_position) = [];
% FDRpass_2(:,p_position) = [];
% 
% FDRpass_3(:,p_position) = [];
% FDRpass_3(:,p_position) = [];
% 
% FDRpass_4(:,p_position) = [];
% FDRpass_4(:,p_position) = [];











% nWrite = ['qpass_T1.xlsx'];
% 
% xlswrite(nWrite,FDRpass_1);
% 
% nWrite = ['qpass_T2.xlsx'];
% 
% xlswrite(nWrite,FDRpass_2);
% 
% nWrite = ['qpass_T3.xlsx'];
% 
% xlswrite(nWrite,FDRpass_3);
% 
% nWrite = ['qpass_T4.xlsx'];
% 
% xlswrite(nWrite,FDRpass_4);

% nWrite = ['qpass_ALL.xlsx'];
% 
% xlswrite(nWrite,FDRpass_all);



%Now, calculate the means and SD of gloged data


result4 = readmatrix('glog_predictedresult.csv');



Col = size (result4,2);
Row = size (result4,1);
result4 = num2cell(result4);



Means_SDs = zeros(n+1,21);
t = 1;
for t = 1:n
    Means_SDs(t,1) = result4{t,1};
end


for s = 1:n
    QCs = zeros(QC,1);
    m = 1;
    for m = 1:QC
        QCs(m,1) = result4{s,m+1};
    end
    QCs = QCs(~isnan(QCs));
    average = mean(QCs);
    Means_SDs(s,2) = average;
    S = std(QCs);
    Means_SDs(s,3) = S;
end

for s = 1:n
    E33s = zeros(E33,1);
    m = 1;
    for m = 1:E33
        E33s(m,1) = result4{s,m+QC+1};
    end
    E33s = E33s(~isnan(E33s));
    average = mean(E33s);
    Means_SDs(s,4) = average;
    S = std(E33s);
    Means_SDs(s,5) = S;
end


for s = 1:n
    E32s = zeros(E32,1);
    m = 1;
    for m = 1:E32
        E32s(m,1) = result4{s,m+QC+E33+1};
    end
    E32s = E32s(~isnan(E32s));
    average = mean(E32s);
    Means_SDs(s,6) = average;
    S = std(E32s);
    Means_SDs(s,7) = S;
end


for s = 1:n
    E31s = zeros(E31,1);
    m = 1;
    for m = 1:E31
        E31s(m,1) = result4{s,m+QC+E33+E32+1};
    end
    E31s = E31s(~isnan(E31s));
    average = mean(E31s);
    Means_SDs(s,8) = average;
    S = std(E31s);
    Means_SDs(s,9) = S;
end


for s =1:n
    E23s = zeros(E23,1);
    m = 1;
    for m = 1:E23
        E23s(m,1) = result4{s,m+QC+E33+E32+E31+1};
    end
    E23s = E23s(~isnan(E23s));
    average = mean(E23s);
    Means_SDs(s,10) = average;
    S = std(E23s);
    Means_SDs(s,11) = S;
end

for s =1:n
    E22s = zeros(E22,1);
    m = 1;
    for m = 1:E22
        E22s(m,1) = result4{s,m+QC+E33+E32+E31+E23+1};
    end
    E22s = E22s(~isnan(E22s));
    average = mean(E22s);
    Means_SDs(s,12) = average;
    S = std(E22s);
    Means_SDs(s,13) = S;
end

for s =1:n
    E21s = zeros(E21,1);
    m = 1;
    for m = 1:E21
        E21s(m,1) = result4{s,m+QC+E33+E32+E31+E23+E22+1};
    end
    E21s = E21s(~isnan(E21s));
    average = mean(E21s);
    Means_SDs(s,14) = average;
    S = std(E21s);
    Means_SDs(s,15) = S;
end

for s =1:n
    E13s = zeros(E13,1);
    m = 1;
    for m = 1:E13
        E13s(m,1) = result4{s,m+QC+E33+E32+E31+E23+E22+E21+1};
    end
    E13s = E13s(~isnan(E13s));
    average = mean(E13s);
    Means_SDs(s,16) = average;
    S = std(E13s);
    Means_SDs(s,17) = S;
end


for s =1:n
    E12s = zeros(E12,1);
    m = 1;
    for m = 1:E12
        E12s(m,1) = result4{s,m+QC+E33+E32+E31+E23+E22+E21+E13+1};
    end
    E12s = E12s(~isnan(E12s));
    average = mean(E12s);
    Means_SDs(s,18) = average;
    S = std(E12s);
    Means_SDs(s,19) = S;
end

for s =1:n
    E11s = zeros(E11,1);
    m = 1;
    for m = 1:E11
        E11s(m,1) = result4{s,m+QC+E33+E32+E31+E23+E22+E21+E13+E12+1};
    end
    E11s = E11s(~isnan(E11s));
    average = mean(E11s);
    Means_SDs(s,20) = average;
    S = std(E11s);
    Means_SDs(s,21) = S;
end

for s =1:n
    CT3s = zeros(CT3,1);
    m = 1;
    for m = 1:CT3
        CT3s(m,1) = result4{s,m+QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+1};
    end
    CT3s = CT3s(~isnan(CT3s));
    average = mean(CT3s);
    Means_SDs(s,22) = average;
    S = std(CT3s);
    Means_SDs(s,23) = S;
end

for s =1:n
    CT2s = zeros(CT2,1);
    m = 1;
    for m = 1:CT2
        CT2s(m,1) = result4{s,m+QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+1};
    end
    CT2s = CT2s(~isnan(CT2s));
    average = mean(CT2s);
    Means_SDs(s,24) = average;
    S = std(CT2s);
    Means_SDs(s,25) = S;
end



for s =1:n
    CT1s = zeros(CT1,1);
    m = 1;
    for m = 1:CT1
        CT1s(m,1) = result4{s,m+QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+CT2+1};
    end
    CT1s = CT1s(~isnan(CT1s));
    average = mean(CT1s);
    Means_SDs(s,26) = average;
    S = std(CT1s);
    Means_SDs(s,27) = S;
end


for s =1:n
    CT0s = zeros(CT0,1);
    m = 1;
    for m = 1:CT0
        CT0s(m,1) = result4{s,m+QC+E33+E32+E31+E23+E22+E21+E13+E12+E11+CT3+CT2+CT1+1};
    end
    CT0s = CT0s(~isnan(CT0s));
    average = mean(CT0s);
    Means_SDs(s,28) = average;
    S = std(CT0s);
    Means_SDs(s,29) = S;
end


s = 1;
for s = 1:n
    
    m = 2;
    for m = 2:29
        if isnan(Means_SDs (s,m))
            Means_SDs (s,m) = 0;
    end
    end
end


Means_SDs(all(Means_SDs==0,2),:) = [];
Means_SDs(:,all(Means_SDs==0,1))= [];

nWrite = ['means&SD_afterGLOG.xlsx'];

xlswrite(nWrite,Means_SDs);




%Combine the significant peaks information with gloged data
File = dir(fullfile(Path,'*xlsx'));
FileNames = {File.name};

result5 = xlsread('glog_predictedresult.csv');
FileCount = size(File,1);
for fc = 1: FileCount
    Filename = FileNames(fc);
    Fname = cell2mat(Filename);
    
    if Fname(end-10:end) == "_qpass.xlsx"
    SigPeak_Info = xlsread(Fname,'Sheet1');
    
    SigPeak_Number = size(SigPeak_Info,1);
    Info = zeros(SigPeak_Number,29);
 for s = 1:SigPeak_Number
     mass = SigPeak_Info(s,1);
     m = 1;
     for m = 1:n
     if (Means_SDs(m,1)-mass) <= 0.000002
         Info(s,:) = Means_SDs(m,:); 

        end
     end 
 end
 
 for s = 1:SigPeak_Number
      mass = SigPeak_Info(s,1);
     m = 1;
     for m = 1:n
     if (result5(m,1)-mass) <= 0.000002
         Sigpeaksinformation(s,:) = result5(m,:); 

        end
     end 
 end
    Fname(end-9:end) = [];
     nWrite = [Fname 'sigpeaks_afterGLOG.xlsx'];

    xlswrite(nWrite,Info); % Use Info for means and SD, use sigpeaksinformation for individual information
    
    
    SampleType = size(Info,2);

    NumSigPeak = size(Info,1);
    
    for t = 1:NumSigPeak

    for m = 1:SampleType-1
    data(t,m) = Info(t,m+1);
    end
    end
    realdata = zeros(NumSigPeak,(SampleType-1)/2);
    for t = 1:(SampleType-1)/2    
    for m = 1:NumSigPeak
    realdata(m,t) = data(m,2*t-1);
    end
    end
    for t = 1:NumSigPeak
    m = 1;
    Zscore_intensity = zeros(1,(SampleType-1)/2);
    for m = 1:(SampleType-1)/2
    Zscore_intensity(1,m)  = realdata(t,m);
    end
    Zscore_rows = zscore(Zscore_intensity);
    Zscore_final(t,:) = Zscore_rows(1,:);
    end
    
    nWrite = [Fname 'sigpeaks_zscored.xlsx'];

    xlswrite(nWrite,Zscore_final);
    
    Sigpeaksinformation = [];
    Zscore_final = [];
    data = [];
    
    
    Info = [];
    end
end


File = dir(fullfile(Path,'*xlsx'));
FileNames = {File.name};

FileCount = size(File,1);
for fc = 1: FileCount
    Filename = FileNames(fc);
    Fname = cell2mat(Filename);
    
    if Fname(7:end) == "_sigpeaks_afterGLOG.xlsx"
    findID = xlsread(Fname,'Sheet1');
    Fname_ID = Fname(1:7);
    Fname_ID = append(Fname_ID,'qvalues');
    IDs = readcell(Fname_ID );
    
    SigPeak_Number = size(findID,1);
    
    IDobtained = num2cell(zeros(SigPeak_Number,1));
 for s = 1:SigPeak_Number
     mass = findID(s,1);
     m = 1;
     for m = 1:n
         mID = cell2mat(IDs(m,1));
     if (mID-mass) <= 0.000002
         IDobtained(s) = IDs(m,5); 

     end
     end 
 end
 
   Fname(end) = [];
   Fname(end) = [];
   Fname(end) = [];
   Fname(end) = [];
   
    nWrite = [Fname '_IDs.xlsx'];

    xlswrite(nWrite,IDobtained);
    
    
    
    
    end
end


% 
% SigPeak_Info = xlsread('qpass_ALL.xlsx','Sheet1');%Change the significant peaks document here
% SigPeak_Number = size(SigPeak_Info,1);
% Info = zeros(SigPeak_Number,21);
% s = 1;
% 
%  for s = 1:SigPeak_Number
%      mass = SigPeak_Info(s,1);
%      m = 1;
%      for m = 1:n
%      if (Means_SDs(m,1)-mass) <= 0.000002
%          Info(s,:) = Means_SDs(m,:); 
% 
%         end
%      end 
%  end
%  result5 = xlsread('glog_predictedresult.csv');
%  for s = 1:SigPeak_Number
%       mass = SigPeak_Info(s,1);
%      m = 1;
%      for m = 1:n
%      if (result5(m,1)-mass) <= 0.000002
%          Sigpeaksinformation(s,:) = result5(m,:); 
% 
%         end
%      end 
%  end
%  
%  nWrite = ['sigpeaks_afterGLOG.xlsx'];
% 
% xlswrite(nWrite,Sigpeaksinformation);


 
%  
%  
% % %Heat maps 
% SampleType = size(Info,2);
% 
% NumSigPeak = size(Info,1);
% yvalues = ["Quality Controls","Exposed TP4","Exposed TP3","Exposed TP2","Exposed TP1","Control TP4","Control TP3","Control TP2","Control TP1","Time Point 0"];
% xvalues = cell(1,NumSigPeak);
% t =1 ;
% for t = 1:NumSigPeak
% xvalues{1,t} = num2str(Info(t,1));
% end
% data = zeros(NumSigPeak,20);
% t = 1;
% for t = 1:NumSigPeak
% m = 1;
% for m = 1:20
% data(t,m) = Info(t,m+1);
% end
% end
% t= 1;
% for t = 1:NumSigPeak
% m = 1;
% for m = 1:20
% data(t,m) = Info(t,m+1);
% end
% end
% Sample_Type = yvalues;
% Metabolite_MZ = xvalues;
% realdata = zeros(NumSigPeak,10);
% t = 1;
% for t = 1:10
% m = 1;
% for m = 1:NumSigPeak
% realdata(m,t) = data(m,2*t-1);
% end
% end
% realdata = realdata';
% 
% %below is the part for putting endpoints in order
% Value_Order = realdata';
% realdata = realdata';
% t =1 ;
% for t = 1:NumSigPeak
% m = 1;
% row_intensity = zeros (1,10);
% for m = 1:10
% row_intensity(1,m) = realdata(t,m);
% end
% [Xsorted,Xidx] = sort(row_intensity);
% [Xsorted,Xidx2] = sort(Xidx);
% Value_Order (t,:) = Xidx2(1,:);
% end
% Value_Order = Value_Order';
% Zscore_final = zeros(NumSigPeak,10);
% 
% t =1 ;
% for t = 1:NumSigPeak
% m = 1;
% Zscore_intensity = zeros(1,10);
% for m = 1:10
% Zscore_intensity(1,m)  = realdata(t,m);
% end
% Zscore_rows = zscore(Zscore_intensity);
% Zscore_final(t,:) = Zscore_rows(1,:);
% end
% ZZZ = Zscore_final';
% 
% %clustergram
% cgo = clustergram (ZZZ);
% Metabolites = Info(:,1);
% Metabolites = Metabolites';
% set(cgo,'ColumnLabels',Metabolites);
% ytest = cellstr(yvalues);
% set(cgo,'RowLabels',ytest);
% set(cgo,'Linkage','complete','Dendrogram',3)
% set(cgo,'Colormap',autumn);
% 
% metabolitesIDs = Sigpeaksinformation(:,1);
% metabolitesIDs = metabolitesIDs';
% ZZZs = [metabolitesIDs;ZZZ];
% 
% nWrite = ['Zscored_sigdata.xlsx'];
% 
% xlswrite(nWrite,ZZZs);
% plot(cgo);
% saveas(gcf,'filename1.png') 
