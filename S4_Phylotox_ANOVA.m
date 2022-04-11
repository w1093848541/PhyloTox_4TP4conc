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

AnovaTPpass = 1;
AnovaECpass = 1;
AnovaIApass = 1;


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
    
    A_Data = [E33t;E32t;E31t;E23t;E22t;E21t;E13t;E12t;E11t;CT3t;CT2t;CT1t];
    ValidE33 = length(E33t);
    ValidE32 = length(E32t);
    ValidE31 = length(E31t);
    ValidE23 = length(E23t);
    ValidE22 = length(E22t);
    ValidE21 = length(E21t);
    ValidE13 = length(E13t);
    ValidE12 = length(E12t);
    ValidE11 = length(E11t);
    
    ValidCT3 = length(CT3t);
    ValidCT2 = length(CT2t);
    ValidCT1 = length(CT1t);
    
    ValidE3 = ValidE33 +ValidE32 + ValidE31;
    ValidE2 = ValidE23 +ValidE22 + ValidE21;
    ValidE1 = ValidE13 +ValidE12 + ValidE11;
    ValidCT = ValidCT3 +ValidCT2 + ValidCT1;
    
    
    AnovaTP = zeros(ValidE3+ValidE2+ValidE1+ValidCT,1);
    AnovaEC = zeros(ValidE3+ValidE2+ValidE1+ValidCT,1);
    
    for TP = 1:ValidE33
        AnovaTP(TP,1) = 3;
    end
    for TP = ValidE33+1:ValidE33+ValidE32
        AnovaTP(TP,1) = 2;
    end
    for TP = ValidE33+ValidE32+1:ValidE33+ValidE32+ValidE31
        AnovaTP(TP,1) = 1;
    end
    for TP = ValidE33+ValidE32+ValidE31+1:ValidE33+ValidE32+ValidE31+ValidE23
        AnovaTP(TP,1)= 3;
    end
    for TP = ValidE33+ValidE32+ValidE31+ValidE23+1:ValidE33+ValidE32+ValidE31+ValidE23+ValidE22
        AnovaTP(TP,1) = 2;
    end
    for TP = ValidE33+ValidE32+ValidE31+ValidE23+ValidE22+1:ValidE33+ValidE32+ValidE31+ValidE23+ValidE22+ValidE21
        AnovaTP(TP,1) = 1;
    end
    for TP = ValidE33+ValidE32+ValidE31+ValidE23+ValidE22+ValidE21+1:ValidE33+ValidE32+ValidE31+ValidE23+ValidE22+ValidE21+ValidE13
        AnovaTP(TP,1) = 3;
    end
    for TP = ValidE33+ValidE32+ValidE31+ValidE23+ValidE22+ValidE21+ValidE13+1:ValidE33+ValidE32+ValidE31+ValidE23+ValidE22+ValidE21+ValidE13+ValidE12
        AnovaTP(TP,1) = 2;
    end
    for TP = ValidE33+ValidE32+ValidE31+ValidE23+ValidE22+ValidE21+ValidE13+ValidE12+1:ValidE33+ValidE32+ValidE31+ValidE23+ValidE22+ValidE21+ValidE13+ValidE12+ValidE11
        AnovaTP(TP,1) = 1;
    end
    for TP = ValidE33+ValidE32+ValidE31+ValidE23+ValidE22+ValidE21+ValidE13+ValidE12+ValidE11+1:ValidE33+ValidE32+ValidE31+ValidE23+ValidE22+ValidE21+ValidE13+ValidE12+ValidE11+ValidCT3
        AnovaTP(TP,1) = 3;
    end   
    for TP = ValidE33+ValidE32+ValidE31+ValidE23+ValidE22+ValidE21+ValidE13+ValidE12+ValidE11+ValidCT3+1:ValidE33+ValidE32+ValidE31+ValidE23+ValidE22+ValidE21+ValidE13+ValidE12+ValidE11+ValidCT3+ValidCT2
        AnovaTP(TP,1) = 2;
    end       
    for TP = ValidE33+ValidE32+ValidE31+ValidE23+ValidE22+ValidE21+ValidE13+ValidE12+ValidE11+ValidCT3+ValidCT2+1:ValidE33+ValidE32+ValidE31+ValidE23+ValidE22+ValidE21+ValidE13+ValidE12+ValidE11+ValidCT3+ValidCT2+ValidCT1
        AnovaTP(TP,1) = 1;
    end          
   
    
    
    for EC = 1:ValidE3
        AnovaEC(EC,1) = 4;
    end
    for EC = 1+ValidE3:ValidE3+ValidE2
        AnovaEC(EC,1) = 3;
    end
    for EC = 1+ValidE3+ValidE2:ValidE3+ValidE2+ValidE1
        AnovaEC(EC,1) = 2;
    end   
    for EC = 1+ValidE3+ValidE2+ValidE1:ValidE3+ValidE2+ValidE1+ValidCT
        AnovaEC(EC,1) = 1;
    end  
    
    
    
     AnovaP = anovan(A_Data,{AnovaTP AnovaEC},'model',2,'varnames',{'Time Point','Exposed vs. Control'},'display','off');
    
    if AnovaP(1,1) <= 0.05
        ts = 1;
        for ts = 1:s
            AnovaTPresult(AnovaTPpass,ts) = cell2mat(result3(t,ts));
            AnovaTPresult(AnovaTPpass,s+1) = AnovaP(1,1);
        end
        AnovaTPpass = AnovaTPpass+1;
    end
    
    if AnovaP(2,1) <= 0.05
        ts = 1;
        for ts = 1:s
            AnovaECresult(AnovaECpass,ts) = cell2mat(result3(t,ts));
            AnovaECresult(AnovaECpass,s+1) = AnovaP(2,1);
        end
        AnovaECpass = AnovaECpass+1;
    end
       
    if AnovaP(3,1) <= 0.05
        ts = 1;
        for ts = 1:s
            AnovaIAresult(AnovaIApass,ts) = cell2mat(result3(t,ts));
            AnovaIAresult(AnovaIApass,s+1) = AnovaP(3,1);
        end
        AnovaIApass = AnovaIApass+1;
    end
end

nWrite = ['TP_ANOVA.xlsx'];

xlswrite(nWrite,AnovaTPresult);

nWrite = ['EC_ANOVA.xlsx'];

xlswrite(nWrite,AnovaECresult);

nWrite = ['IA_ANOVA.xlsx'];

xlswrite(nWrite,AnovaIAresult);

Anova_all = [AnovaTPresult;AnovaECresult;AnovaIAresult];
p_position = size(Anova_all,2);
Anova_all(:,p_position) = [];

nWrite = ['all_ANOVA.xlsx'];
xlswrite(nWrite,Anova_all);

result4 = readmatrix('glog_predictedresult.csv');



Col = size (result4,2);
Row = size (result4,1);
result4 = num2cell(result4);



% Means_SDs = zeros(n+1,21);
% t = 1;
% for t = 1:n
%     Means_SDs(t,1) = result4{t,1};
% end
% 
% 
% for s = 1:n
%     QCs = zeros(QC,1);
%     m = 1;
%     for m = 1:QC
%         QCs(m,1) = result4{s,m+1};
%     end
%     QCs = QCs(~isnan(QCs));
%     average = mean(QCs);
%     Means_SDs(s,2) = average;
%     S = std(QCs);
%     Means_SDs(s,3) = S;
% end
% 
% for s = 1:n
%     E4s = zeros(E4,1);
%     m = 1;
%     for m = 1:E4
%         E4s(m,1) = result4{s,m+QC+1};
%     end
%     E4s = E4s(~isnan(E4s));
%     average = mean(E4s);
%     Means_SDs(s,4) = average;
%     S = std(E4s);
%     Means_SDs(s,5) = S;
% end
% 
% 
% for s = 1:n
%     E3s = zeros(E3,1);
%     m = 1;
%     for m = 1:E3
%         E3s(m,1) = result4{s,m+QC+E4+1};
%     end
%     E3s = E3s(~isnan(E3s));
%     average = mean(E3s);
%     Means_SDs(s,6) = average;
%     S = std(E3s);
%     Means_SDs(s,7) = S;
% end
% 
% 
% for s = 1:n
%     E2s = zeros(E2,1);
%     m = 1;
%     for m = 1:E2
%         E2s(m,1) = result4{s,m+QC+E4+E3+1};
%     end
%     E2s = E2s(~isnan(E2s));
%     average = mean(E2s);
%     Means_SDs(s,8) = average;
%     S = std(E2s);
%     Means_SDs(s,9) = S;
% end
% 
% 
% for s =1:n
%     E1s = zeros(E1,1);
%     m = 1;
%     for m = 1:E1
%         E1s(m,1) = result4{s,m+QC+E4+E3+E2+1};
%     end
%     E1s = E1s(~isnan(E1s));
%     average = mean(E1s);
%     Means_SDs(s,10) = average;
%     S = std(E1s);
%     Means_SDs(s,11) = S;
% end
% 
% 
% for s = 1:n
%     C4s = zeros(C4,1);
%     m = 1;
%     for m = 1:C4
%         C4s(m,1) = result4{s,m+QC+E4+E3+E2+E1+1};
%     end
%     C4s = C4s(~isnan(C4s));
%     average = mean(C4s);
%     Means_SDs(s,12) = average;
%     S = std(C4s);
%     Means_SDs(s,13) = S;
% end
% 
% 
% for s = 1:n
%     C3s = zeros(C3,1);
%     m = 1;
%     for m = 1:C3
%         C3s(m,1) = result4{s,m+QC+E4+E3+E2+E1+C4+1};
%     end
%     C3s = C3s(~isnan(C3s));
%     average = mean(C3s);
%     Means_SDs(s,14) = average;
%     S = std(C3s);
%     Means_SDs(s,15) = S;
% end
% 
% 
% for s = 1:n
%     C2s = zeros(C2,1);
%     m = 1;
%     for m = 1:C2
%         C2s(m,1) = result4{s,m+QC+E4+E3+E2+E1+C4+C3+1};
%     end
%     C2s = C2s(~isnan(C2s));
%     average = mean(C2s);
%     Means_SDs(s,16) = average;
%     S = std(C2s);
%     Means_SDs(s,17) = S;
% end
% 
% 
% for s = 1:n
%     C1s = zeros(C1,1);
%     m = 1;
%     for m = 1:C1
%         C1s(m,1) = result4{s,m+QC+E4+E3+E2+E1+C4+C3+C2+1};
%     end
%     C1s = C1s(~isnan(C1s));
%     average = mean(C1s);
%     Means_SDs(s,18) = average;
%     S = std(C1s);
%     Means_SDs(s,19) = S;
% end
% 
% 
% for s = 1:n
%     C0s = zeros(C0,1);
%     m = 1;
%     for m = 1:C0
%         C0s(m,1) = result4{s,m+QC+E4+E3+E2+E1+C4+C3+C2+C1+1};
%     end
%     C0s = C0s(~isnan(C0s));
%     average = mean(C0s);
%     Means_SDs(s,20) = average;
%     S = std(C0s);
%     Means_SDs(s,21) = S;
% end
% 
% s = 1;
% for s = 1:n
%     
%     m = 2;
%     for m = 2:21
%         if isnan(Means_SDs (s,m))
%             Means_SDs (s,m) = 0;
%     end
%     end
% end
% 
% 
% Means_SDs(all(Means_SDs==0,2),:) = [];
% Means_SDs(:,all(Means_SDs==0,1))= [];
% 
% nWrite = ['means&SD_afterGLOG.xlsx'];
% 
% xlswrite(nWrite,Means_SDs);
       

 result5 = xlsread('glog_predictedresult.csv');
 Means_SDs  = xlsread('means&SD_afterGLOG.xlsx');

SigPeak_Info = xlsread('TP_ANOVA.xlsx','Sheet1');%Change the significant peaks document here
SigPeak_Number = size(SigPeak_Info,1);
Info = zeros(SigPeak_Number,29);
s = 1;

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
         SigpeaksinformationTP(s,:) = result5(m,:); 

        end
     end 
 end
 
 nWrite = ['TP_ANOVA_sigpeaks_afterGLOG.xlsx'];

xlswrite(nWrite,Info);

 
 
 
 
%Heat maps 
% SampleType = size(Info,2);
% 
% NumSigPeak = size(Info,1);
% yvalues = ["QCs","Exposed TP4","Exposed TP3","Exposed TP2","Exposed TP1","Control TP4","Control TP3","Control TP2","Control TP1","TP0"];
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


SigPeak_Info = xlsread('EC_ANOVA.xlsx','Sheet1');%Change the significant peaks document here
SigPeak_Number = size(SigPeak_Info,1);
Info = zeros(SigPeak_Number,29);
s = 1;

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
         SigpeaksinformationEC(s,:) = result5(m,:); 

        end
     end 
 end
 
 nWrite = ['EC_ANOVA_sigpeaks_afterGLOG.xlsx'];

xlswrite(nWrite,Info);

 
 
 
 
% %Heat maps 
% SampleType = size(Info,2);
% 
% NumSigPeak = size(Info,1);
% yvalues = ["QCs","Exposed TP4","Exposed TP3","Exposed TP2","Exposed TP1","Control TP4","Control TP3","Control TP2","Control TP1","TP0"];
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


SigPeak_Info = xlsread('IA_ANOVA.xlsx','Sheet1');%Change the significant peaks document here
SigPeak_Number = size(SigPeak_Info,1);
Info = zeros(SigPeak_Number,29);
s = 1;

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
         SigpeaksinformationIA(s,:) = result5(m,:); 

        end
     end 
 end
 
 nWrite = ['IA_ANOVA_sigpeaks_afterGLOG.xlsx'];

xlswrite(nWrite,Info);






File = dir(fullfile(Path,'*xlsx'));
FileNames = {File.name};

result5 = xlsread('glog_predictedresult.csv');
FileCount = size(File,1);
for fc = 1: FileCount
    Filename = FileNames(fc);
    Fname = cell2mat(Filename);
    
    if Fname(3:end) == "_ANOVA_sigpeaks_afterGLOG.xlsx"
    Info = xlsread(Fname,'Sheet1');
    
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
    Fname(9:end) = [];
    nWrite = [Fname '_sigpeaks_zscored.xlsx'];

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
    
    if Fname(3:end) == "_ANOVA_sigpeaks_afterGLOG.xlsx"
    findID = xlsread(Fname,'Sheet1');
    
    IDs = readcell('T3_CV1_qvalues.xlsx' );%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
 
 
   Fname(end-5:end) = [];
    
    nWrite = [Fname '_IDs.xlsx'];

    xlswrite(nWrite,IDobtained);
    
    
    
    
    end
end






 
 
 
 
% %Heat maps 
% SampleType = size(Info,2);
% 
% NumSigPeak = size(Info,1);
% yvalues = ["QCs","Exposed TP4","Exposed TP3","Exposed TP2","Exposed TP1","Control TP4","Control TP3","Control TP2","Control TP1","TP0"];
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
%     
    