function [pdf,c]=blmlestimate(X,xTest,fb,options)
% function to calculate blml density estimate using blmlTrivial algorithm
% For details see Agarwal R, Chen Z, Sarma SV, A Novel Nonparametric Maximum Likelihood
% Estimator for Probability Density Functions. IEEE TPAMI 2016.
%Inputs- X: observations nxndim
%        xTest: test points where pdf is evaluated mxndim
%        fb: 1xndim vector of cutoff frequencies
%        options: option structure containing options.memSize, MaxIter, and Tolerance
%Outputs- pdf: Joint PDF values evaluated at points in xgrid.
%
%Examples-
%         2-d pdf:
%         [xgrid{1:2}]=ndgrid(-5:0.05:5);
%         pdf=blmlestimate([randn(4000,1) randn(4000,1)],[xgrid{1}(:) xgrid{2}(:)],[1 1],options);
%         surf(xgrid{1},xgrid{2},reshape(pdf,size(xgrid{1})))
% author: Rahul Agarwal. rahul.jhu@gmail.com

[n,ndim]=size(X);

if n<options.memSize
    [t1s, t2s]=meshgrid(1:n);
    s=ones(size(t1s));
    for i=1:ndim
        x=X(:,i);
        s=fb(i)*s.*sinc(fb(i)*(x(t2s)-x(t1s)));

    end
    c0=ones(n,1);
    options=optimset('Jacobian','on','TolX',options.TolX,'TolFun',options.TolFun,'MaxIter',options.MaxIter,'Display','off');
    c=fsolve(@(c)fc(c,s),c0/n,options);
    k=c'*s*c;
    c=c/sqrt(k);
else
    c=solve2f(X,fb,options)/length(X);
end

pdf=zeros(size(xTest,1),1);


for i=1:n

rx=ones(size(xTest,1),1);
    for j=1:ndim
    rx=fb(j)*rx.*sinc(fb(j)*(xTest(:,j)-X(i,j)));
    end
pdf=pdf+c(i)*rx;
end
pdf=pdf.^2;

