disp('See http://www.gaussianprocess.org/gpml/code/matlab/doc/ for details.')
disp('Hit any key to continue...'); pause

disp(' '); disp('clear all, close all')
clear all, close all
write_fig = 0;
disp(' ')

disp('n1 = 80; n2 = 40;                   % number of data points from each class')
n1 = 80; n2 = 40;
disp('S1 = eye(2); S2 = [1 0.95; 0.95 1];           % the two covariance matrices')
S1 = eye(2); S2 = [1 0.95; 0.95 1];
disp('m1 = [0.75; 0]; m2 = [-0.75; 0];                            % the two means')
m1 = [0.75; 0]; m2 = [-0.75; 0];
disp(' ')

disp('x1 = bsxfun(@plus, chol(S1)''*gpml_randn(0.2, 2, n1), m1);')
x1 = bsxfun(@plus, chol(S1)'*gpml_randn(0.2, 2, n1), m1);
disp('x2 = bsxfun(@plus, chol(S2)''*gpml_randn(0.3, 2, n2), m2);')         
x2 = bsxfun(@plus, chol(S2)'*gpml_randn(0.3, 2, n2), m2);         
disp(' ')

disp('n1 = 800; n2 = 400;                   % number of data points from each class')
tn1 = 100; tn2 = 100;
%disp('S1 = eye(2); S2 = [1 0.95; 0.95 1];           % the two covariance matrices')
tS1 = eye(2); tS2 = [1 0.95; 0.95 1];
disp('m1 = [0.75; 0]; m2 = [-0.75; 0];                            % the two means')
tm1 = [0.75; 0]; tm2 = [-0.75; 0];
disp(' ')

disp('x1 = bsxfun(@plus, chol(S1)''*gpml_randn(0.2, 2, n1), m1);')
tx1 = bsxfun(@plus, chol(tS1)'*gpml_randn(0.2, 2, tn1), tm1);
disp('x2 = bsxfun(@plus, chol(S2)''*gpml_randn(0.3, 2, n2), m2);')         
tx2 = bsxfun(@plus, chol(tS2)'*gpml_randn(0.3, 2, tn2), tm2);         
disp(' ')

tx = [tx1 tx2]'; ty = [-ones(1,tn1) ones(1,tn2)]';


disp('x = [x1 x2]''; y = [-ones(1,n1) ones(1,n2)]'';')
x = [x1 x2]'; y = [-ones(1,n1) ones(1,n2)]';
figure(6)
disp('plot(x1(1,:), x1(2,:), ''b+''); hold on;');
plot(x1(1,:), x1(2,:), 'b+', 'MarkerSize', 12); hold on
disp('plot(x2(1,:), x2(2,:), ''r+'');');
plot(x2(1,:), x2(2,:), 'r+', 'MarkerSize', 12);
disp(' ')

disp('[t1 t2] = meshgrid(-4:0.1:4,-4:0.1:4);')
[t1 t2] = meshgrid(-4:0.1:4,-4:0.1:4);
disp('t = [t1(:) t2(:)]; n = length(t);               % these are the test inputs')
t = [t1(:) t2(:)]; n = length(t);               % these are the test inputs
disp('tmm = bsxfun(@minus, t, m1'');')
tmm = bsxfun(@minus, t, m1');
disp('p1 = n1*exp(-sum(tmm*inv(S1).*tmm/2,2))/sqrt(det(S1));')
p1 = n1*exp(-sum(tmm*inv(S1).*tmm/2,2))/sqrt(det(S1));
disp('tmm = bsxfun(@minus, t, m2'');')
tmm = bsxfun(@minus, t, m2');
disp('p2 = n2*exp(-sum(tmm*inv(S2).*tmm/2,2))/sqrt(det(S2));')
p2 = n2*exp(-sum(tmm*inv(S2).*tmm/2,2))/sqrt(det(S2));
set(gca, 'FontSize', 24)
disp('contour(t1, t2, reshape(p2./(p1+p2), size(t1)), [0.1:0.1:0.9])')
contour(t1, t2, reshape(p2./(p1+p2), size(t1)), [0.1:0.1:0.9])
[c h] = contour(t1, t2, reshape(p2./(p1+p2), size(t1)), [0.5 0.5]);
set(h, 'LineWidth', 2)
colorbar
grid
axis([-4 4 -4 4])
if write_fig, print -depsc f6.eps; end
disp(' '); disp('Hit any key to continue...'); pause

disp(' ')
disp('meanfunc = @meanConst; hyp.mean = 0;')
meanfunc = @meanConst; hyp.mean = 0;
disp('covfunc = @covSEard;   hyp.cov = log([1 1 1]);')
%covfunc = @covSEard;   hyp.cov = log([1 1 1]);
covfunc = @covLIN;   hyp.cov = [];
disp('likfunc = @likErf;')
likfunc = @likErf;
disp(' ')

disp('hyp = minimize(hyp, @gp, -40, @infEP, meanfunc, covfunc, likfunc, x, y);')
hyp = minimize(hyp, @gp, -40, @infEP, meanfunc, covfunc, likfunc, x, y);
disp('[a b c d lp] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, x, y, t, ones(n, 1));')
[a b c d lp] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, x, y, t, ones(n,1));
disp(' ')
figure(7)
set(gca, 'FontSize', 24)
disp('plot(x1(1,:), x1(2,:), ''b+''); hold on')
plot(x1(1,:), x1(2,:), 'b+', 'MarkerSize', 12); hold on
disp('plot(x2(1,:), x2(2,:), ''r+'')')
plot(x2(1,:), x2(2,:), 'r+', 'MarkerSize', 12)
disp('contour(t1, t2, reshape(exp(lp), size(t1)), 0.1:0.1:0.9)')
contour(t1, t2, reshape(exp(lp), size(t1)), 0.1:0.1:0.9)
[c h] = contour(t1, t2, reshape(exp(lp), size(t1)), [0.5 0.5]);
set(h, 'LineWidth', 2)
colorbar
grid
axis([-4 4 -4 4])
if write_fig, print -depsc f7.eps; end
disp(' '); disp('Hit any key to continue...'); pause

disp('large scale classification using the FITC approximation')
disp('[u1,u2] = meshgrid(linspace(-2,2,5)); u = [u1(:),u2(:)];')
[u1,u2] = meshgrid(linspace(-2,2,5)); u = [u1(:),u2(:)]; clear u1; clear u2
disp('nu = size(u,1);')
nu = size(u,1);
disp('covfuncF = {@covFITC, {covfunc}, u};')
covfuncF = {@covFITC, {covfunc}, u};
disp('inffunc = @infFITC_Laplace;')
inffunc = @infFITC_EP;                     % one could also use @infFITC_Laplace
disp('hyp = minimize(hyp, @gp, -40, inffunc, meanfunc, covfuncF, likfunc, x, y);')
hyp = minimize(hyp, @gp, -40, inffunc, meanfunc, covfuncF, likfunc, x, y);
disp('[a b c d lp] = gp(hyp, inffunc, meanfunc, covfuncF, likfunc, x, y, t, ones(n,1));')
%[a b c d lp] = gp(hyp, inffunc, meanfunc, covfuncF, likfunc, x, y, t, ones(n,1));
[a b c d lp] = gp(hyp, inffunc, meanfunc, covfuncF, likfunc, x, y, tx, ones(n,1));
disp(' ')
figure(8)
set(gca, 'FontSize', 24)
disp('plot(x1(1,:), x1(2,:), ''b+''); hold on')
plot(x1(1,:), x1(2,:), 'b+', 'MarkerSize', 12); hold on
disp('plot(x2(1,:), x2(2,:), ''r+'')')
plot(x2(1,:), x2(2,:), 'r+', 'MarkerSize', 12)
disp('contour(t1, t2, reshape(exp(lp), size(t1)), 0.1:0.1:0.9)')
%contour(tx1, tx2, reshape(exp(lp), size(tx1)), 0.1:0.1:0.9)
%[c h] = contour(tx1, tx2, reshape(exp(lp), size(tx1)), [0.5 0.5]);
%set(h, 'LineWidth', 2)
%disp('plot(u(:,1),u(:,2),''ko'', ''MarkerSize'', 12)')
%plot(u(:,1),u(:,2),'ko', 'MarkerSize', 12)
%colorbar
%grid
%axis([-4 4 -4 4])
%if write_fig, print -depsc f8.eps; end
%disp(' ')

correctPredictions = 0;
for m = 1:200
  if m>=1 && m<101
      if a(m)<0
          correctPredictions = correctPredictions + 1;
          
      end
      
  elseif m>=101
      if a(m)>0
          correctPredictions = correctPredictions + 1;
      end
  end
end

accuracy = correctPredictions/200;

accuracy

class1 = abs((a<0) - (ty<0)) == 0 ; % if 1 its correct
class2 = abs((a<0) - (ty<0)) ~= 0; % if 1 its incorrect

correctPredictions = lp(class1);
wrongPredictions = lp(class2);

pCorrect = exp(correctPredictions);
pWrong = exp(wrongPredictions);

for i=1:size(correctPredictions)
   cEntr(i) = -1*(pCorrect(i)*log2(pCorrect(i)) + (1-pCorrect(i))*log2(1-pCorrect(i))); 
end

for i=1:size(wrongPredictions)
    wEntr(i) = -1*((pWrong(i)*log2(pWrong(i))) + ((1-pWrong(i))*log2(1-pWrong(i))));
end


cEntr = cEntr';
wEntr = wEntr';

figure(9)
hist(cEntr, 10);
figure(10)
hist(wEntr, 10);


%templTree = templateTree('MaxNumSplits',15,'Surrogate','on');
ClassTreeEns = fitensemble(x, y, 'Bag', 200, 'Tree','Type', 'Classification');
[labelsPredict score] = predict(ClassTreeEns, tx);
classificationRate = sum((labelsPredict == ty)) / size(ty, 1)

classC = (labelsPredict == ty);
classW = (labelsPredict ~= ty);

t=1;
for p=1:200
   if classC(p) == 1 
    correctBoostingProbability(t) = max(score(p,:));
    t = t + 1;
   end
end
k=1;
for p=1:200
   if classW(p) == 1 
    wrongBoostingProbability(k) = max(score(p,:));
    k = k + 1;
   end
end

correctBoostingProbability = correctBoostingProbability';
wrongBoostingProbability = wrongBoostingProbability';

for i=1:size(correctBoostingProbability)
   correctBoostingEntr(i) = -1*(correctBoostingProbability(i)*log2(correctBoostingProbability(i)) + (1-correctBoostingProbability(i))*log2(1-correctBoostingProbability(i))); 
end

for i=1:size(wrongBoostingProbability)
    wrongBoostingEntr(i) = -1*((wrongBoostingProbability(i)*log2(wrongBoostingProbability(i))) + ((1-wrongBoostingProbability(i))*log2(1-wrongBoostingProbability(i))));
end


correctBoostingEntr = correctBoostingEntr';
wrongBoostingEntr = wrongBoostingEntr';

figure(12)
hist(correctBoostingEntr, 10);
figure(13)
hist(wrongBoostingEntr, 10);




