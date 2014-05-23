function reg = MultiPolyRegressV2(Data,R,PW,PV)
%   Multivariable polynomial regression analysis. Output 'reg' is a struct with
%   the fitted constants, coressponding terms and R-Square. 
%
%   reg = MultiPolyRegress(Data,R,PW) performs multivariable polymonial 
%   regression analysis on row stacked dimensional data matrix Data. Data is 
%   an m-by-n matrix where m and n are the number of dimensions and the 
%   number of data points. R is the n-by-1 response vector and PW is the degree
%   of the polynomial fit. 
%   
%   reg = MultiPolyRegress(Data,R,PW,PV) restricts individual dimensions of
%   Data to particular powers PV in the polynomial expansion. PV is an
%   m-by-1 vector.
%   
%   Uses the MATLAB Statistical and Symbolic Toolboxes.
%
%   Author : Ahmet Cecen

    % Align Data
    if size(Data,2)>size(Data,1)
        Data=Data';
    end
    
    % DegreeWorks
    if nargin == 3
        PV = repmat(PW,[1,size(Data,2)]);
    end
    
    % Function Parameters
    NData = size(Data,1);
    NVars = size(Data,2);
    RowMultiB = '1';
    RowMultiC = '1';
    Lim = max(PV);
    
    % Initialize
    A=zeros(Lim^NVars,NVars);

    % Create Colums Corresponding to Mathematical Base
    for i=1:NVars
        A(:,i)=mod(floor((1:Lim^NVars)/Lim^(i-1)),Lim);
    end

    % Flip - Reduce - Augment
    A=fliplr(A); A=A(sum(A,2)<=Lim,:); Ab=diag(repmat(Lim,[1,NVars])); A=[A;Ab];

    % Degree Conditionals
    for i=1:NVars
        A=A(A(:,i)<=PV(i),:);
    end

    % Build Framework
    B=sym(zeros(size(A,1),NVars));
    for i=1:NVars
        B(:,i)=sym(['x',num2str(i)]);
        RowMultiB=strcat(RowMultiB,['.*B(:,',num2str(i),')']);
        RowMultiC=strcat(RowMultiC,['.*C(:,',num2str(i),')']);
    end
    
    % Create a Legend for Coefficient Correspondence
    B=B.^A; Legend = eval(RowMultiB); %#ok<NASGU>
    
    % Allocate
    NLegend = length(Legend);
    Scores = zeros(NData,NLegend);
    
    % Compose
    for i=1:NData
        current=repmat(Data(i,:),[NLegend,1]);
        C=current.^A; %#ok<NASGU>
        Scores(i,:) = eval(RowMultiC);
    end

    % Regress
    [b,bint,r,rint,stats] = regress(R,Scores,0.05);
    RSq=stats(1);
    
    reg = struct('PowerMatrix',A,'Scores',Scores,'Coefficients', b,'Residuals', r, 'Legend', Legend, 'RSquare', RSq,'RInt',rint,'BInt',bint);
end
    
    