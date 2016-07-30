function pdf=kdeestimate(X,xTest,fc)

%Inputs: X: observations nxndim
%        xTest: ndimx1 cell array containing points where the pdf is
%        evaluated. Elements of the cell array are either the combination
%        of m-d matrices or scalars. Where 0<=m<=ndim.
%        fc: constant in sigma=0.4/fc(j)/n^0.2 where sigma is std of the
%        kernel used
%Outputs: pdf: PDF values evaluated at points in xgrid.
%Examples-
%         2-d pdf:
%         [xgrid{1:2}]=ndgrid(-5:0.05:5);
%         pdf=kdeestimate([randn(1000,1) randn(1000,1)],[xgrid{1}(:) xgrid{2}(:)],[1 1]);
%         surf(xgrid{1},xgrid{2},reshape(pdf,size(xgrid{1})))
% author: Rahul Agarwal. rahul.jhu@gmail.com

[n,ndim]=size(X);
pdf=zeros(size(xTest,1),1);

for i=1:n
rx=ones(size(xTest,1),1);
    for j=1:ndim
    rx=rx.*normpdf((xTest(:,j)-X(i,j)),0,0.4/fc(j)/n^(0.2));
    end
    
pdf=pdf+rx;
end

pdf=pdf/n;