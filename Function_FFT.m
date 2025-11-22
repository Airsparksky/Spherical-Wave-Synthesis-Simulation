function  [f0,fp,fn]=Function_FFT(Value,m) 
n=m/2;
f1=fft(Value,m);
f0=f1(1);
for i=1:n
    fp(i)=f1(i+1);%取值从第一个开始到中间的值，得到f1的前一半
    fn(i)=f1(m+1-i);%取值从最后一个开始到中间的值，得到f1的后一半
end

delphi=2*pi/m;
f01=sum(Value)*delphi/(2*pi);
sp1=0;
sp2=0;
sp3=0;
sp0=0;
sn1=0;
sn2=0;
sn3=0;
for nl=0:(length(Value)-1)
    sp0=sp0+Value(nl+1);
    sp1=sp1+(Value(nl+1))*(cos(nl*delphi)-j*sin(nl*delphi));
    sp2=sp2+(Value(nl+1))*(cos(2*nl*delphi)-j*sin(nl*2*delphi));
    sp3=sp3+(Value(nl+1))*(cos(3*nl*delphi)-j*sin(nl*3*delphi));
    
    sn1=sn1+(Value(nl+1))*(cos(nl*delphi)+j*sin(nl*delphi));
    sn2=sn2+(Value(nl+1))*(cos(2*nl*delphi)+j*sin(nl*2*delphi));
    sn3=sn3+(Value(nl+1))*(cos(3*nl*delphi)+j*sin(nl*3*delphi));
end
sp11=sp1*delphi/(2*pi);
sp21=sp2*delphi/(2*pi);
sp31=sp3*delphi/(2*pi);

sn11=sn1*delphi/(2*pi);
sn21=sn2*delphi/(2*pi);
sn31=sn3*delphi/(2*pi);
f0=f0/(m);
fp=fp/(m);
fn=fn/(m);