clear all;
basemva=100;
baseKV=11;
 nbus=4;
 %        lp  lq  r       x      ysh   tap
linedata=[1  2   2       -8      0    1
          1  3   1       -4      0    1
          2  3   0.666   -2.664  0    1
          2  4   1       -4      0    1
          3  4   2       -8      0    1];
nline=length(linedata(:,1))
j=sqrt(-1);
i=sqrt(1);

for k=1:nline
    lp(k)=linedata(k,1);
    lq(k)=linedata(k,2);
    r(k)=linedata(k,3);
    x(k)=linedata(k,4);
    ysh(k)=linedata(k,5);
    a(k)=linedata(k,6);
    z(k)=r(k)+j*x(k);
    y(k)=1/z(k);
end
lp
lq
r
x
ysh
a

y
 
ybus=zeros(nbus,nbus);
yln=zeros(nbus,nbus);

for k=1:nline
    ylp(k)= [1/(a(k)^2)-1/(a(k))]*y(k);
    ylq(k)= [1-1/a(k)]*y(k);
    yt(k)=y(k)/a(k);
end
ylp
ylq
y
 for k=1:nline
ybus(lp(k),lq(k)) = ybus(lp(k),lq(k)) - yt(k);
ybus(lq(k),lp(k)) = ybus(lp(k),lq(k));
ybus(lp(k),lp(k)) = ybus(lp(k),lp(k)) + yt(k) + ylp(k) + j*ysh(k);
ybus(lq(k),lq(k)) = ybus(lq(k),lq(k)) + yt(k) + ylq(k) + j*ysh(k);


 end 
 ybus