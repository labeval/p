clc;
clear all;
%First element is added to reference  
%both from and to bus nodes should not be new nodes
%New node is written under From and old node in To
%     element  From   To    Z
% data=[ 1        1     0    0.25;
%        2        2     1    0.15;
%        3        3     1    0.15;
%        4        2     0    0.25;
%        5        2     3    0.15 ];
%     element  From   To    Z 
data=[ 1        1     0    0.6;
       2        2     1    0.25;
       3        3     2    0.5;
       4        3     0    0.5 ];
   
element=data(:,1);   
nbr=length(data(:,1));                   %No. of Branches
from=data(:,2);                          %From Bus
to=data(:,3);                            %To Bus
Zb=data(:,4);                            %Impedence     
n=max(max(from),max(to));                %No. of Busses
Zbus=zeros(n,n);
for i=1:nbr
%Modification-1
%A New bus is added to reference bus
if(element(i)==1)
Zbus=Zb(i);
continue
end
%Modification-2
%A New bus is added to old bus other than reference bus
if(from(i)~=0 && to(i)~=0)
if (from(i)>to(i))
k=to(i);
new=from(i);
for j=1:2
Zbus(j,new)=Zbus(j,k);
Zbus(new,j)=Zbus(k,j);
end
Zbus(new,new)=Zbus(k,k)+Zb(i);
continue
end
end
%Modification-3
%A existing bus is added to reference bus
if(to(i)==0)
old=from(i);
m1=Zbus(old,old)+Zb(i);
ztemp=(1/m1)*Zbus(:,old)*Zbus(old,:);
Zbus=Zbus-ztemp;
continue
end
%Modification-4
%An existing bus is added to another existing bus other than reference bus
if(from(i)~=0 && to(i)~=0)
a=from(i);
b=to(i);
m2=Zb(i)+Zbus(a,a)+Zbus(b,b)-(2*Zbus(a,b));
ztemp=(1/m2)*((Zbus(:,a)-(Zbus(:,b)))*((Zbus(a,:))-(Zbus(b,:))));
Zbus=Zbus-ztemp;
continue
end      
end
fprintf('Z-Bus\n');
disp(Zbus);