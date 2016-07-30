function pdf=blmlestimatequick(X,xTest,fb,options)
% function to calculate blml density estimate using blmlquick algorithm
%  For details see Agarwal R, Chen Z, Sarma SV, A Novel Nonparametric Maximum Likelihood
%  Estimator for Probability Density Functions. IEEE TPAMI 2016.
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


[nall,ndim]=size(X);
res=1/((nall)^0.25)./fb;
X=floor(X./(ones(size(X,1),1)*res)+1/2).*(ones(size(X,1),1)*res);
[x,~,ic]=unique(X,'rows');
n=length(x);
xn=zeros(n,1);
for i=1:n
    xn(i)=length(find(ic==i));
end

if n<options.memSize
    [t1s,t2s]=meshgrid(1:n);
    s=ones(size(t1s));
    for i=1:ndim
        xi=x(:,i);
        s=fb(i)*s.*sinc(fb(i)*(xi(t2s)-xi(t1s)));
    end
    s=s*diag(xn);
    c0=ones(n,1);
    options = optimset('Jacobian','on','TolX',options.TolX,'TolFun',options.TolFun,'MaxIter',options.MaxIter,'Display','off');
    c = fsolve(@(c)fc(c,s),c0/n,options);
    c=c/sqrt(nall);
else
     c=solve2fquick(x,fb,xn,options)/nall;
end

pdf=zeros(size(xTest,1),1);

for i=1:n

rx=ones(size(xTest,1),1);
    for j=1:ndim
    rx=fb(j)*rx.*sinc(fb(j)*(xTest(:,j)-x(i,j)));
    end
    
pdf=pdf+c(i)*xn(i)*rx;
end
pdf=pdf.^2;



