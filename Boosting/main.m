% Change the filenames if you've saved the files under different names
% On some platforms, the files might be saved as 
% train-images.idx3-ubyte / train-labels.idx1-ubyte
trainImages = loadMNISTImages('train-images-idx3-ubyte');
trainLabels = loadMNISTLabels('train-labels-idx1-ubyte');

testImages = loadMNISTImages('t10k-images-idx3-ubyte');
testLabels = loadMNISTLabels('t10k-labels-idx1-ubyte');

realTrainImages = transpose(trainImages);
realTestImages = transpose(testImages);

Ensemble = fitensemble(realTrainImages(1:1000, :), trainLabels(1:1000, :), 'AdaBoostM2', 100, 'Tree');
[results,prob] = predict(Ensemble, realTestImages);

rsLoss = resubLoss(Ensemble, 'Mode','Cumulative');
plot(rsLoss);
xlabel('Number of learning cycles');
ylabel('Resub loss');

correctPredictions = 0;
for m = 1:10000
  if testLabels(m) == results(m)
      correctPredictions = correctPredictions + 1;
  end
end

accuracy = correctPredictions / 10000;

accuracy