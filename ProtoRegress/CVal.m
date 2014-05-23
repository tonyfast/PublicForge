function [CV, CVCoeff, CVLeg]=CVal(PCA,R,NPC,PW)
% Cross-Validation wrapper for polynomial regression. Uses symbolic toolbox to create the legend, might need to make a legend-less version.
% The inputs are the inputs to MultiPolyRegress, CV is the Cross Validated score vector, CVCoeff stores the coefficients at each iteration
% (comparison can be used as some measure of stability), CVLeg is legend at each iteration, for correspondence between coefficients.
n=size(PCA,1);
PC=PCA(:,1:NPC);
for i=1:n
    current=circshift(PC,[-i 0]);
    R1=circshift(R,-i);
    Data=current(1:(n-1),:);
    R2=R1(1:(n-1));
    reg = MultiPolyRegressV2(Data,R2,PW);
    PM=reg.PowerMatrix;
    current=repmat(PC(i,:),[size(PM,1),1]);
    C=current.^PM;
    nn=size(C,2);
    cc1=ones(size(C,1),1);
    for k=1:nn
        cc1=cc1.*C(:,k);
    end
    CV(i)=cc1'*reg.Coefficients;
    CVCoeff{i}=reg.Coefficients;
    CVLeg{i}=reg.Legend;
end