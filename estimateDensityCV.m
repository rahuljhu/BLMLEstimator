function [pdf,fb]=estimateDensityCV(x,xTest,type,options,fb)
%Inputs- x: i.i.d observations nxndim
%        xTest: test points where pdf is evaluated mxndim
%        type (optional): 'blml', 'blmlq' 'kde' or 'hist' 
%               For details see Agarwal R, Chen Z, Sarma SV, A Novel Nonparametric Maximum Likelihood
%               Estimator for Probability Density Functions. IEEE TPAMI 2016.  
%        options (optional): option structure containing options.memSize(size of side of max square matrix memory can hold, default 4000,
%                             if n, nbins(for blmlq)>options.memSize, solve2f (based on gradient descent) is used to solve rho_n(c)=0), 
%                            MaxIter (default 1000), TolX (default 1e-9), TolFun (default 1e-9) 
%        fb (optional): 1xndim vector of cutoff frequencies
%
%Outputs- pdf: Joint PDF values evaluated at points in xTest.
%          fb: Estimated cut-off frequencies 1xndim
% 
%Examples-
%         2-d pdf:
%         [xgrid{1:2}]=ndgrid(-5:0.05:5);
%         [pdf,fb]=estimateDensityCV(randn(1000,2),[xgrid{1}(:) xgrid{2}(:)]);
%         surf(xgrid{1},xgrid{2},reshape(pdf,size(xgrid{1})))
% author: Rahul Agarwal. rahul.jhu@gmail.com

if nargin <=2
    type='blmlq';
end

if nargin<=3
    options.MaxIter=1000;
    options.TolX=1e-9;
    options.TolFun=1e-9;
    options.memSize=4000;
end

if nargin<=4
    xTrain=x(1:floor(end*0.5),:);
    xCV=x(floor(end*0.5)+1:end,:);
    fbv=(2.^(-4:0.1:1000))';%1./(0.5:0.1:4);
    fbv=fbv*(1./(quantile(x,0.75)-quantile(x,0.25)));
    logl=-Inf(length(fbv),1);
    for j=1:length(fbv)
        if (strcmp(type,'blmlq'))
            logl(j)=mean(log(blmlestimatequick(xTrain,xCV,fbv(j,:),options)));
        elseif (strcmp(type,'blml'))
            logl(j)=mean(log(blmlestimate(xTrain,xCV,fbv(j,:),options)));
        elseif (strcmp(type,'hist'))
            logl(j)=mean(log(histEstimate(xTrain,xCV,fbv(j,:))));
        else
            logl(j)=mean(log(kdeestimate(xTrain,xCV,fbv(j,:))));
        end
        [~,ind]=max(logl);
        if ind+10<j
            break;
        end
    end
    [~,j]=max(logl);
    fb=fbv(j,:);
end



if(strcmp(type,'blmlq'))
    pdf=blmlestimatequick(x,xTest,fb,options);
elseif(strcmp(type,'blml'))
    pdf=blmlestimate(x,xTest,fb,options);
elseif (strcmp(type,'hist'))
    pdf=histEstimate(x,xTest,fb);
else
    pdf=kdeestimate(x,xTest,fb);
end
