function pdf=histEstimate(X,xTest,fb)
% function for estimating histogram
%Inputs- X: observations nxndim
%        xTest: test points where pdf is evaluated mxndim
%        fb: 1xndim vector of cutoff frequencies
%Outputs- pdf: Joint PDF values evaluated at points in xgrid.
%
%Examples-
%         2-d pdf:
%         [xgrid{1:2}]=ndgrid(-5:0.05:5);
%         pdf=histEstimate([randn(4000,1) randn(4000,1)],[xgrid{1}(:) xgrid{2}(:)],[1 1]);
%         surf(xgrid{1},xgrid{2},reshape(pdf,size(xgrid{1})))
% author: Rahul Agarwal. rahul.jhu@gmail.com

res=1./fb;
X=floor(X./(ones(size(X,1),1)*res)+1/2).*(ones(size(X,1),1)*res);
[x,~,ic]=unique(X,'rows');
n=size(x,1);
xn=zeros(n,1);
for i=1:n
    xn(i)=length(find(ic==i));
end
pn=xn/size(X,1)/prod(res);
xTest=floor(xTest./(ones(size(xTest,1),1)*res)+1/2).*(ones(size(xTest,1),1)*res);
[Xgrid,~,xgic]=unique(xTest,'rows');
pdf=zeros(size(Xgrid,1),1);
for i=1:size(Xgrid,1)
    ind=ones(n,1);
    for j=1:length(res)
        ind=ind==1 & x(:,j)==Xgrid(i,j);
    end
    p=pn(ind);
    if ~isempty(p)
        pdf(i)=p;
    else
        pdf(i)=0;
    end
end
pdf=pdf(xgic);
