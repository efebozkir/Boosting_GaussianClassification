% EFE BOZKIR, TU Munich
% Reading the MNIST Dataset
trainImages = loadMNISTImages('train-images-idx3-ubyte');
trainLabels = loadMNISTLabels('train-labels-idx1-ubyte');

testImages = loadMNISTImages('t10k-images-idx3-ubyte');
testLabels = loadMNISTLabels('t10k-labels-idx1-ubyte');

realTrainImages = transpose(trainImages);
realTestImages = transpose(testImages);

% Training operation
Ensemble = fitensemble(realTrainImages(1:1000, :), trainLabels(1:1000, :), 'AdaBoostM2', 100, 'Tree');
%Ensemble = fitensemble(realTrainImages(1:60000, :), trainLabels(1:60000, :), 'Bag', 1000, 'Tree','Type', 'Classification');
% Testing the test set
[results,prob] = predict(Ensemble, realTestImages);

% Plot the resubstritution loss with learning cycles
rsLoss = resubLoss(Ensemble, 'Mode','Cumulative');
plot(rsLoss);
xlabel('Number of learning cycles');
ylabel('Resub loss');

% Calculation of accuracy of the system
correctPredictions = 0;
for m = 1:10000
  if testLabels(m) == results(m)
      correctPredictions = correctPredictions + 1;
  end
end

accuracy = correctPredictions / 10000;

accuracy
