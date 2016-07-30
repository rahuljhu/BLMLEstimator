% script for showing example run of estimateDensityCV
% author: Rahul Agarwal. rahul.jhu@gmail.com
%% 1-d example
clear
err=zeros(100,3);
dx=0.001;
xgrid=(-5:dx:5)';
pdftrue=normpdf(xgrid,0,1);
for i=1:100
x=randn(1000,1);
[pdfK,fbK]=estimateDensityCV(x,xgrid,'kde');
err(i,1)=sum((pdfK-pdftrue).^2)*dx;
[pdf,fb]=estimateDensityCV(x,xgrid,'blmlq');
err(i,2)=sum((pdf-pdftrue).^2)*dx;
[pdf,fbH]=estimateDensityCV(x,xgrid,'hist');
err(i,3)=sum((pdf-pdftrue).^2)*dx;
display(i);
end
mean(err)
[~,pK]=ttest(err(:,1),err(:,2))
[~,pH]=ttest(err(:,3),err(:,2))

%% 2-d example
clear
err=zeros(100,3);
dx=0.05;
[xgrid{1:2}]=ndgrid((-5:dx:5)');
pdftrue=mvnpdf([xgrid{1}(:) xgrid{2}(:)]);
for i=1:10
x=randn(4000,2);
[pdfK,fbK]=estimateDensityCV(x,[xgrid{1}(:) xgrid{2}(:)],'kde');
err(i,1)=sum((pdfK-pdftrue).^2)*dx^2;
[pdf,fb]=estimateDensityCV(x,[xgrid{1}(:) xgrid{2}(:)],'blmlq');
err(i,2)=sum((pdf-pdftrue).^2)*dx^2;
[pdf,fbH]=estimateDensityCV(x,[xgrid{1}(:) xgrid{2}(:)],'hist');
err(i,3)=sum((pdf-pdftrue).^2)*dx^2;
display(i);
end
mean(err)
[~,pK]=ttest(err(:,1),err(:,2))
[~,pH]=ttest(err(:,3),err(:,2))