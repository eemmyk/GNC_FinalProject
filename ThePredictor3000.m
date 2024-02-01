close all;
clear;

load traningData_Earth.csv

dataSize = size(traningData_Earth);

% Data structure: 
% 1          2   3    4   5       6        7   8    9   10      11           14           13     
% a_initial, P1, Tp1, e1, omega1, a_final, P2, Tp2, e2, omega2, dateOptimal, tof_optimal, deltaV_date


rng(1);

outPutInds = [11, 12, 13];

normallyNormalize = [1, 2, 6, 7];
piNormalize = [5, 10];
noNormalize = [4, 9];
relNormalize = [3, 8];

outputCount = numel(outPutInds);

inputCount = numel([normallyNormalize, piNormalize, noNormalize, relNormalize]);

inputs = zeros(dataSize(1), inputCount);
outputs = zeros(dataSize(1), outputCount); 

%Relative Normalization for 3 and 8
inputs(:, 3) = traningData_Earth(:, 3) ./ traningData_Earth(:, 2);
inputs(:, 8) = traningData_Earth(:, 8) ./ traningData_Earth(:, 7);

%Normalize inputs
for i = normallyNormalize
    min_input = min(traningData_Earth(:, i));
    span_input = max(traningData_Earth(:, i)) - min(traningData_Earth(:, i));

    inputs(:, i) = (traningData_Earth(:, i) - min_input) / span_input;
end

for i = piNormalize
    inputs(:, i) = traningData_Earth(:, i) / (2*pi);
end

for i = noNormalize
    inputs(:, i) = traningData_Earth(:, i);
end

%Process outputs

outputs(:, 1) = traningData_Earth(:, outPutInds(1));
outputs(:, 2) = log(1+traningData_Earth(:, outPutInds(1)));
outputs(:, 3) = traningData_Earth(:, outPutInds(1));

y_in_Date = outputs(:, 1)';
y_in_TOF = outputs(:, 2)';
y_in_DV = outputs(:, 3)';

%Transpose inputs and outputs

x_in = inputs';
y_in = outputs(:, 1)';

maxLayerSize1 = 50;
minLayerSize1 = 10;

maxLayerSize2 = 20;
minLayerSize2 = 20;

vec_rmse_train = zeros(outputCount, maxLayerSize1, maxLayerSize2);
vec_rmse_val = zeros(outputCount, maxLayerSize1, maxLayerSize2);

figure;
hold on;

legendTexts = {};

for j = minLayerSize2:maxLayerSize2
    for i = minLayerSize1:maxLayerSize1
        %Define the ANN
        hiddenLayerSize = [i, j];
        
        net = fitnet(hiddenLayerSize);
        net.divideParam.trainRatio = 70/100;
        net.divideParam.valRatio = 30/100;
        net.divideParam.testRatio = 0/100;
        [net, tr] = train(net, x_in, y_in);
            
        yTrain = net(x_in(:, tr.trainInd));
        yTrainTrue = y_in(:, tr.trainInd);
        
        vec_rmse_train(1, i, j) = sqrt(mean((yTrain(1, :) - yTrainTrue(1, :)).^2));
%         vec_rmse_train(2, i, j) = sqrt(mean((yTrain(2, :) - yTrainTrue(2, :)).^2));
%         vec_rmse_train(3, i, j) = sqrt(mean((yTrain(3, :) - yTrainTrue(3, :)).^2));

        yVal = net(x_in(:, tr.valInd));
        yValTrue = y_in(:, tr.valInd);
        
        vec_rmse_val(1, i, j) = sqrt(mean((yVal(1, :) - yValTrue(1, :)).^2));
%         vec_rmse_val(2, i, j) = sqrt(mean((yVal(2, :) - yValTrue(2, :)).^2));
%         vec_rmse_val(3, i, j) = sqrt(mean((yVal(3, :) - yValTrue(3, :)).^2));
    end

    plot(1:maxLayerSize1, vec_rmse_train(1, :, j));
    plot(1:maxLayerSize1, vec_rmse_val(1, :, j));

    drawnow;

    legendTexts(end+1) = {''};
    legendTexts(end+1) = {sprintf("layer2: %.0f", j)};
end

legend(legendTexts);

% 
% yTrain = net(x_in(:, tr.trainInd));
% yTrainTrue = y_in(:, tr.trainInd);
% 
% yVal = net(x_in(:, tr.valInd));
% yValTrue = y_in(:, tr.valInd);
% 
% 
% figure;
% hold on;
% 
% plot(yValTrue, yVal, 'o');
% 
% figure;
% hold on;
% 
% plot(yTrainTrue, yTrain, 'o');

% %Define the ANN
% hiddenLayerSize = [9,2];
% 
% DV_net = fitnet(hiddenLayerSize);
% DV_net.divideParam.trainRatio = 70/100;
% DV_net.divideParam.valRatio = 30/100;
% DV_net.divideParam.testRatio = 0/100;
% [DV_net, tr] = train(DV_net, x_in, y_in_DV);
% 
% yTrain = DV_net(x_in(:, tr.trainInd));
% yTrainTrue = y_in_DV(:, tr.trainInd);
% 
% yVal = DV_net(x_in(:, tr.valInd));
% yValTrue = y_in_DV(:, tr.valInd);
% 
% figure;
% hold on;
% 
% plot(yValTrue, yVal, 'o');
% 
% %% Define final Date model
% 
% %Define the ANN
% hiddenLayerSize = [100,40];
% 
% Date_net = fitnet(hiddenLayerSize);
% Date_net.divideParam.trainRatio = 70/100;
% Date_net.divideParam.valRatio = 30/100;
% Date_net.divideParam.testRatio = 0/100;
% [Date_net, tr] = train(Date_net, x_in, y_in_Date);
% 
% 
% yTrain = Date_net(x_in(:, tr.trainInd));
% yTrainTrue = y_in_Date(:, tr.trainInd);
% 
% yVal = Date_net(x_in(:, tr.valInd));
% yValTrue = y_in_Date(:, tr.valInd);
% 
% 
% figure;
% hold on;
% 
% plot(yValTrue, yVal, 'o');
