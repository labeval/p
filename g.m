clc
clear all

% This is for case-7 and case-8
% MVAb=100;           %Base MVA
% KVb=12.66;          %Base kV

MVAb=100;
KVb=11;
a=1/(MVAb*1000);B=MVAb/KVb^2;

% Define the line data of 5 bus system
          % lp 	lq    G      B    ysh 	tap
linedata =[ 1 	2     2    	-8      0 	1
            1 	3     1 	-4      0	1
            2  	3 	0.666 	-2.664 	0 	1
            2 	4     1 	-4      0 	1
            3 	4     2    	-8      0 	1];

        
P=[0 .5 .4 .3];
Q=[0 .2 .3 .1];

v = [1.06 1 1 1];

from=linedata(:,1);                          %From Bus
to=linedata(:,2);                            %To Bus
nbus=max(max(from),max(to))  ;            %No. of Busses
        

nline=length(linedata(:,1));
j=sqrt(-1);
i=sqrt(1); 

for k=1:nline
lp(k)=linedata(k,1);
lq(k)=linedata(k,2);
r(k)=linedata(k,3);     % original data may be given in Y or Z form
x(k)=linedata(k,4);     % original data may be given in Y or Z form
ysh(k)=linedata(k,5);
a(k)=linedata(k,6);
y(k)=(r(k)+j*x(k));     % y(k) is taken becoz data is given in Y form
%y(k)=1/z(k);
end

ybus=zeros(nbus,nbus);
yln =zeros(nbus,nbus);

% PI METHOD FOR OFF-NOMINAL ADMITTANCE OF TRANSFORMER
for k=1:nline
ylp(k)=[1/(a(k)^2)-1/a(k)]*y(k);
ylq(k)=[1-1/a(k)]*y(k);
y(k)=y(k)/a(k); 
end 

for k=1:nline
ybus(lp(k),lq(k))=ybus(lp(k),lq(k))-y(k); % This for off diagonal value
ybus(lq(k),lp(k))=ybus(lp(k),lq(k));
ybus(lp(k),lp(k))=ybus(lp(k),lp(k))+y(k)+ylp(k)+j*ysh(k);
ybus(lq(k),lq(k))=ybus(lq(k),lq(k))+y(k)+ylq(k)+j*ysh(k);
end

ybus;


alpha=1.6;
Er=0.0001;
%total_bus=4;

for Iteration=1:1
    for p=1:nbus
        if p==1
            V(p)=v(1);  
        else
            X1= (1/ybus(p,p))*(((-P(p)+Q(p)*1i)/conj(v(p))));
            
            X2=0;
        for q=1:p-1       
            X2=X2+(1/ybus(p,p))*(ybus(p,q)*v(q));
        end  
            X3=0;
        for q=p+1:nbus        
            X3=X3+(1/ybus(p,p))*((ybus(p,q)*v(q)));
        end
            V(p)=X1-X2-X3;
        end
        p;
        V;
        v;
             
        %for gauss-seidel method
        Vacc(p)=v(p)+alpha*(V(p)-v(p));
       %error
        ER(p)=sqrt((abs(Vacc(p))-abs(v(p)))^2);  
  
        %v(p)=Vacc(p); %remove % for GSM and put % for GM
        %v(1)=1.06;
    end
    ER;
    if max(ER)<Er
         break
     end
    %v(p)=Vacc(p);
    Iteration
    %for gauss method
    v= V;
    v(1)=1.06;
    Vx=V ;   
end
Vx
for z=1:nbus
Voltage(z)=abs(Vx(z));
end
Voltage


%V4=(1/(3-12*1j))*((-0.3+0.1*1j)-((-1+4*1j)*(1.01899-0.046208*1j))-((-2+8*1j)*(0.99059-0.0467968*1j))) 

            