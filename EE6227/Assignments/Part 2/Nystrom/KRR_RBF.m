function [acc,Trainningtime,Testingtime] = KRR_RBF(TrainX,TrainY,TestX,TestY,option)

kernel = option.kernel;
kernel_param = option.kernel_param;
lambda = option.lambda;

Ulabel = unique(TrainY);
Nsample = size(TrainX,1);
m = round(sqrt(Nsample)); % sqrt of num of samples (default)

%% train
Trainningtime = 0;

tic
trainY_temp = oneVrestCoding(TrainY,Ulabel);

TDindx = 1:size(TrainX,1);
Indx = randsample(TDindx,m);
trainX_sam = TrainX(Indx,:);

K_b = kernel_matrix(TrainX,kernel,kernel_param,trainX_sam);
K_hat = kernel_matrix(trainX_sam,kernel,kernel_param);
K_hat = pinv(K_hat);

inver_matrix = (eye(Nsample,Nsample)-K_b*pinv(lambda*eye(m,m)+K_hat*(K_b'*K_b))*K_hat*K_b')/lambda;
c = inver_matrix*trainY_temp;

Trainingtime_temp = toc;
Trainningtime = Trainningtime + Trainingtime_temp;

%% test
Ntest = size(TestX,1);
Testingtime = 0;

tic
Kt = kernel_matrix(TrainX,kernel,kernel_param,TestX);
Yt = Kt'*c;

Y = oneVrestDecoding(Yt,Ulabel);

Testingtime_temp = toc;
Testingtime = Testingtime + Testingtime_temp;


acc = length(find(Y==TestY))/Ntest;

end

