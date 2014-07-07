M=100;
l=fzero(@(x)2^x-2*M-1, M); 
N=2^(ceil(l));
cVector=zeros; mm=[0:1:2*M];mmInv=[2*M:1:1]
cVector=[exp(1i*pi.*mm.^2*theta) zeros(1,N-4*M-1) exp(1i*pi.*mmInv.^2*theta) ];

