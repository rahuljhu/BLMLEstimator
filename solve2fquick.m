function [c,ndlogl]=solve2fquick(X,fb,xn,c0)
%Inputs- X: observations nxndim
%        fb: ndimx1 vector of cutoff frequencies
%        c0: nx1 vector of labels +-1 indicating class label of each data point
%Outputs- c: solution of rho(c)=0.
%author: Rahul Agarwal. rahul.jhu@gmail.com
[n,ndim]=size(X);
nall=sum(xn);

if nargin==3
    c0=ones(n,1);
else
    c0=c0;
end

a0=zeros(n,1);
for i=1:n
    s=xn;
    for j=1:ndim
        x=X(:,j);
        s=fb(j)*s.*sinc(fb(j)*(x-x(i)));
    end
    s(i)=0;
    a0(i)=s'*c0;
    
end
c1=(3*c0+(sign(c0).*sqrt(a0.^2+4*prod(fb)*nall*xn)-a0)/2/prod(fb)./xn)/4;

while sum((c1-c0).^2)>1e-18
    c0=c1;
    a0=zeros(n,1); 
    for i=1:n
        s=xn;
        for j=1:ndim
            x=X(:,j);
            s=fb(j)*s.*sinc(fb(j)*(x-x(i)));
        end
        s(i)=0;
        a0(i)=s'*c0;
    end
    c1=(3*c0+(sign(c0).*sqrt(a0.^2+4*prod(fb)*nall*xn)-a0)/2/prod(fb)./xn)/4;
end

c=c1;
if nargout>1
    ndlogl=0;
for i=1:length(c)
    vec=0;
    s=1;
    for j=1:ndim
        x=X(:,j);
        s=fb(j)*s.*sinc(fb(j)*(x-x(i)));
        vec=vec+ cos(fb(j)*(x-x(i)))./sinc(fb(j)*(x-x(i)))/fb(j);
    end
    vec=s.*vec;
    ndlogl=ndlogl+sum(xn(i)*c(i)*c.*vec.*xn);
end
ndlogl=ndlogl/nall^2;
end
