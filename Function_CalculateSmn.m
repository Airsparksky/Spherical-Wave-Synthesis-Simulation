function  [SD0,Sp,SDp,Sn,SDn]=Function_CalculateSmn(cost,sint,N,M)
%S 是组合式Smn；%S0=0，
%SD 是组合式Smn的导数,SD0对应m=0 ；
%Sp,SDp对应m>0;
%Sn,SDn对应m<0
%cost t对应Theta 
%SD0 ==> n
%Sp,SDp ==> m*n
%Sn,SDn ==>  m*n
if N>=2 && M<=N    
    v=sqrt(3/(8*pi));    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %m=0的情况
    s00=0;
    SD0(1)=-v*sint;   %  先从S'mn的第二项开始，最后S'mn第一项再添加到SD0中  
    Temp_q=0;
    for n=2:N
        qv2=sqrt((2*n+1)/(2*n-1)*(n-1)/(n+1));
        if n==2
            Temp1=s00;
        else
            Temp1=SD0(n-2);
        end        
        SD0(n)=qv2/(n-1)*((2*n-1)*cost*SD0(n-1)-Temp_q*n*Temp1);
        Temp_q=qv2;
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %m=n m>0
    Temp=v/sqrt(2);
    for n=1:N
        Sp(n,n)=Temp;
        SDp(n,n)=n*Temp*cost;
        Temp=Temp*sint*sqrt((2*n+3)/(2*n+2)*n/(n+2));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %m<n m>0
    for m=1:(M-1)
        k=m+1;
        TempSp1=0;
        TempQ=0;
        TempSp2=Sp(m,m);
        for n=k:N
            Q=sqrt((n-1)*(n-m)*(2*n+1)/((n+m)*(2*n-1)*(n+1)));
            Sp(m,n)=Q/(n-m)*((2*n-1)*cost*TempSp2-(n+m-1)*TempQ*TempSp1);
            SDp(m,n)=n*cost*Sp(m,n)-(n+m)*Q*TempSp2;
            TempSp1=TempSp2;
            TempSp2=Sp(m,n);
            TempQ=Q;
        end
    end
else
     disp(' M 应当大于或等于2,且N大于或等于M.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %根据对称性，求相应的m<0的值
TempSp=Sp(1:M,:);
TempSDp=SDp(1:M,:);
if(mod(M,2)) % M ==> odd
   i=1:2:M;
   TempSp(i,:)=-TempSp(i,:);
   TempSDp(i,:)=-TempSDp(i,:);
else  % M ==> even
   i=1:2:M-1;
   TempSp(i,:)=-TempSp(i,:);
   TempSDp(i,:)=-TempSDp(i,:);
end
Sn=TempSp;
SDn=TempSDp;
