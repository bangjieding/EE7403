clc;
clear;

% Load datasets (The string inside the load function is the path to the
% .mat file)

%abalone\bank\car\contrac\hill-valley\led-display\plant-margin\plant-texture 
train_set = load('../MatDataset/hill-valley/hill-valley_Train.mat');
test_set = load('../MatDataset/hill-valley/hill-valley_Test.mat');

% Collect training/testing dataset
trainX = train_set.Data;
trainY = train_set.Label;
testX = test_set.Data;
testY = test_set.Label;

% Parameters
ModelParameters.nTree = 10;                            % Number of trees [int]
ModelParameters.mtry = round(sqrt(size(trainX,2)));     % Number of features to sample [int]


%% Tuning
% Parameters to tune (You can also experiment with different values)
mtry_range = 1:size(trainX,2);
% mtry_range = 13:30;
mtry_range

if size(trainX,2) > 25
    mtry_range = 1:round(sqrt(size(trainX,2))+20);
    mtry_range
end

best_acc = 0;   % Initialisation

% Requried for consistency
s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);

cv_part = cvpartition(trainY,'KFold',4);    % Create indices for training/validation subsets

% Test every configuration
for p1 = 1:numel(mtry_range)
    p1
    TestModelParameters = ModelParameters;
    TestModelParameters.mtry = mtry_range(p1);

    val_acc = zeros(4,1);   % Initialisation
    parfor k = 1:4
        % Collect training/validation sets
        val_trainX = trainX(cv_part.training(k),:);
        val_trainY = trainY(cv_part.training(k),:);
        val_testX = trainX(cv_part.test(k),:);
        val_testY = trainY(cv_part.test(k),:);
        
        % Data Normalisation
        mean_X = mean(val_trainX,1);
        std_X = std(val_trainX);
        std_X(std_X==0) = 1e-4;                 % For numerical stability
        val_trainX = bsxfun(@rdivide,val_trainX-repmat(mean_X,size(val_trainX,1),1),std_X);
        val_testX = bsxfun(@rdivide,val_testX-repmat(mean_X,size(val_testX,1),1),std_X);

        % Training and Testing
        val_acc(k) = RF(val_trainX,val_trainY,val_testX,val_testY,TestModelParameters);

    end

    % Average the validation accuracy
    ValAcc = mean(val_acc);

    % Check if current configuration is the best
    if ValAcc > best_acc
        best_acc = ValAcc;
        best_mtry = mtry_range(p1);
    end
    % parsave(sprintf('output%d.mat', p1), ValAcc, p1)
end

% Use the best settings
ModelParameters.mtry = best_mtry;
ModelParameters.nTree = 100;

%% Evaluation
% Data Normalisation
mean_X = mean(trainX,1);
std_X = std(trainX);
std_X(std_X==0) = 1;                 % For numerical stability
trainX = bsxfun(@rdivide,trainX-repmat(mean_X,size(trainX,1),1),std_X);
testX = bsxfun(@rdivide,testX-repmat(mean_X,size(testX,1),1),std_X);

[TestAcc,TrainingTime,TestingTime,Model] = RF(trainX,trainY,testX,testY,ModelParameters);

function [] = parsave(dir, ValAcc, p1)
    %save x,y in dir
    % so I can save in parfor loop
    save(dir,'ValAcc', 'p1');
end



