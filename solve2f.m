function c=solve2f(X,fb,options)
%Inputs- X: observations nxndim
%        fb: ndimx1 vector of cutoff frequencies
%        options: option structure containing MaxIter, and Tolerance
%Outputs- c: solution of rho(c)=0.
%author: Rahul Agarwal. rahul.jhu@gmail.com
[n,ndim]=size(X);



c0=ones(n,1)/sqrt(n);


a0=zeros(n,1);
for i=1:n
    s=ones(n,1);
    for j=1:ndim
        x=X(:,j);
        s=fb(j)*s.*sinc(fb(j)*(x-x(i)));
    end
    s(i)=0;
    a0(i)=s'*c0;
    
end
c1=(3*c0+(sign(c0).*sqrt(a0.^2+4*prod(fb)*n)-a0)/2/prod(fb))/4;
iter=0;
while sqrt(sum((c1-c0).^2))>options.TolX && iter<options.MaxIter
    c0=c1;
    a0=zeros(n,1); 
    for i=1:n
        s=ones(n,1);
        for j=1:ndim
            x=X(:,j);
            s=fb(j)*s.*sinc(fb(j)*(x-x(i)));
        end
        s(i)=0;
        a0(i)=s'*c0;
    end
    c1=(3*c0+(sign(c0).*sqrt(a0.^2+4*prod(fb)*n)-a0)/2/prod(fb))/4;
    iter=iter+1;
end
if (iter>=options.MaxIter)
    display('Iteration stopped due to MaxIter Limit reached: Increase options.MaxIter')
end
c=c1;
