%% Data pre-Processing 

% ©Chacha Chen

%% ========Part1:  Importing data====================

clear;close all;clc;

test_load;                     %import all data into variable all
import_label_pre;              %import all data into variable ETABM185
all = table2array(all);         
ETABM185 = ETABM185(2:5897,8);  %delete the first row,and extract the 8th (label) column
label = table2cell(ETABM185);  % 5896*1 cell

fprintf('/ndata importing completed!');
%label = unique(ETABM185(:,8));

%% ========Part2: Filtering the datasets==============

data = transpose(all);


index_todele = find(cellfun(@isempty,label));
data(index_todele,:) = [];

label = label(~cellfun(@isempty,label)); % deleting the empty label

%%
 
labels = unique(label);  
index_todele = [];
 
% check if specific label data <10
for i =1:size(labels,1)
     index_t = contains(label,labels(i));      %index of label = labels(i)

     if (nnz(index_t)<10)
         index_todele = [index_todele;find(index_t)];  %deleting those index
     end
     
 end
 
label(index_todele) = [];
data(index_todele,:) = [];
 
labels = unique(label);
 
fprintf('/nDatasets filtering completed!');     


%% PCA on training data

[coeff, score, latent] = pca(data); % coefficients
recon_deg = cumsum(latent)./sum(latent); 
fprintf('PCA completed!');

fprintf(find(recon_deg > 0.8,1));
fprintf(find(recon_deg > 0.9,1));
fprintf(find(recon_deg > 0.95,1));
   
score = score(:,1:1064);
score = zscore(score);
dlmwrite('dataPro.txt',score);

fprintf('z-score completed!');

%resulting a 3k*1k doubel data matrix and a 3k label matrix 98 labels


%%
%model = svmtrain(label_01, data_norm);
%[lbl, acc, dec] = svmpredict(label_01(1:100,:), data_norm(1:100,:), model);
lable_processed = zeros(size(label,1),1);
for i =1:size(labels,1)
    index_t = contains(label,labels(i));
    index_t = find(index_t);
    label_processed(index_t,1) = i;  
 end
%% crossvalidation
% 
%%
% rng(1);
% 
% SVMModel = fitcsvm(data_norm,label_01,'KernelFunction','RBF',...
%     'KernelScale','auto');
% 
% CVSVMModel = crossval(SVMModel);
% 
% classLoss = kfoldLoss(CVSVMModel);

%% extracting all labels about cancer/tumor
% key_word = ["tumor","cancer","adenocarcinoma",...
% "AIDS-KS",...
% "leukemia",...
% "leukaemia",...
% "lymphoma",...
% "Sarcoma",...
% "adenocarcinoma.",...
% "alveolar rhabdomyosarcoma",...
% "carcinoma",...
% "chondroblastoma",...
% "chondromyxoid fibroma",...
% "chordoma",...
% "colon adenocarcinoma",...
% "dedifferentiated chondrosarcoma",...
% "embryonal rhabdomyosarcoma",...
% "fibromatosis",...
% "follicular thyroid adenoma",...
% "follicular thyroid carcinoma",...
% "ganglioneuroblastoma",...
% "ganglioneuroma",...
% "glioblastoma",...
% "grade 2, primary hnscc",...
% "high-stage neuroblastoma",...
% "hlrcc",...
% "iatrogenic-KS, KSHV-",...
% "leiomyosarcoma",...
% "lipoma",...
% "low-stage neuroblastoma",...
% "lung adenocarcinoma",...
% "sarcoma",...
% "myxoid liposarcoma",...
% "neuroblastoma",...
% "neurofibroma",...
% "osteosarcoma",...
% "adenoma",...
% "schwannoma",...
% "t4b",...
% "tendon xanthomas",...
% "uterine fibroid",...
% "well-differentiated liposarcoma","aml",...
% "Classic-KS, HIV-, nodular (late) stage",...
% "Daudi Burkitt's lymphoma"];
% 
% cancer_index=[];
% 
% for i =1:size(key_word,2)
%     tmp = strfind(label,key_word(i));
%     temp = find(~cellfun(@isempty,tmp));
%     cancer_index = [cancer_index;temp];
% end
% 
% cancer_index = unique(cancer_index);  %2831*1 double
% sort(cancer_index);
% 
% label_01 = zeros(size(label));
% 
% for i = 1: size(cancer_index,1)
%     tmp = cancer_index(i);
%     label_01(tmp) = 1;
% end
% 
% fprintf('label generating completed!');