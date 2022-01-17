close all
clc
format loose
%n=24;            % Order of the graph

B= [[0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0]];


N=[2.33333, 2.50000, 3.00000, 2.50000, 3.00000, 2.33333, 3.00000, 3.00000, 3.00000, 3.00000, 2.50000, 2.33333, 2.50000, 2.50000, 2.50000, 2.33333, 2.33333, 2.50000, 2.50000, 2.33333, 2.50000, 2.50000, 2.50000, 2.50000];



n=size(N',1);

A= reshape(B,[n,n]);
if A==A'
      disp('Matrix is Symmetric')
else
      disp('Not Symmetric')
end

l=size(A,1);
e=[];
s=[];
t=[];
k=[];
for i=1:l
    for j=i+1:l
        if A(i,j)==1
            e=[e;[i,j]];
        end
    end
end

m=size(e,1);
j=ones(m,1);

k=sum(A);
k1=k';

q=[];
for i=1:l
    qq=[i k1(i)];
    q=[q;qq];
end

p=[];
for i=1:m
    pp=[e(i,:) k(e(i,1)) k(e(i,2))];
    p=[p;pp];
end

X1=q(:,2);
X2=X1.^3;

A1=p(:,3);
B1=p(:,4);
A2= A1.*B1;
D2= A1.^2 + B1.^2;
E2= sqrt(D2); 
C2= (A1-1).*(B1-1);
B2=A1+B1;


F=sum(X2);                                         %F(G) index


A6=1/sqrt(A2);                                % Ra index
R=sum(A6);

A3=A2.^(1/2);                                  %RR(G) index
RR=sum(A3);
%
A4=B2;                                         %M_{1}(G) index
M1=sum(A4);

A5=A2;                                         %M_{2}(G) index
M2=sum(A5);



%Albertson's Irregularity Index
A7=A1-B1;
A8=abs(A7);
AL=sum(A8);                                    %AL(G)  


%Sigma Index
sigma=F-2*M2;                                 %Sigma(G)


%Degree Variance
variance=M1./n-((2*m)/n)^2;                   %Var(G)


%IRM1 Index
IRM1=sqrt(M1/n)-((2*m)/n);                   %IRM1(G)

%IR1 Index
IR1=sqrt(M1/m)-((2*m)/n);                    %IR1(G)

%IRM2 Index
IRM2=sqrt(M2/n)-((2*m)/n);                    %IRM2(G)

%IR2 Index
IR2=sqrt(M2/m)-((2*m)/n);                    %IR1(G)

%IRC Index
IRC=F-(((2*m)/n)*M1);                        %IRC(G)

%IRF Index
IRF=F-2*M2;                                 %IRF(G)

%IRFW Index
IRFW=F/M2;                                  %IRFW(G)

%IRV Index
IRV=n*M1-(2*m).^2;                          %IRV(G)

%IRA Index
IRA=n-2*R;                                    %IRA(G)

%IRB Index
IRB=M1-2*RR;                                %IRB(G)

%IRL Index
IRL=RR./m-((2*m)/n);                        %IRL(G)
 
%IRDIF Index
A9=A1./B1;
A10=B1./A1;
A11=A9-A10;
A12=abs(A11);
IRDIF=sum(A12);                                  %IRDIF Index


%IRZ Index
A13=log(A1)-log(B1);
A14=abs(A13);
IRZ=sum(A14);                                  %IRZ Index

%IRLU Index
A15=abs(A7);
A16=min(A1,B1);
A17=A15./A16;
IRLU=sum(A17);                                  %IRLU Index

%IRLF Index
A18=sqrt(A2);
A19=A15./A18;
IRLF=sum(A19);                                  %IRLF Index

%IRLA Index
A20=A15./B2;
IRLA=2*sum(A20);                                  %IRLA Index

%IRD1 Index
A21=1+A15;
A22=log(A21);
IRD1=sum(A22);                                  %IRD1 Index

%IRGA Index
A23=2*sqrt(A2);
A24=B2./A23;
A25=log(A24);
IRGA=sum(A25);                                  %IRGA Index


%IRT1 Index
A26=sum(N);
IRT1=((1/n)*A26)-((2*m)/n);                     %IRT1 Index


%IRT2 Index
N1=N';
A27=k1-N1;
A28=abs(A27);
IRT2=sum(A28);                                    %IRT2 Index




disp('Irregularity indices:')
% fprintf('The Ra(G) Index is %4.4f\n',R');
% fprintf('The RR(G) Index is %4.4f\n',RR');
% fprintf('The M_{1}(G) Index is %4.4f\n',M1');
% fprintf('The M_{2}(G) Index is %4.4f\n',M2');
% fprintf('The F(G) Index is %4.4f\n',F');
%fprintf('The AL(G) Index is %4.4f\n',AL');
fprintf('The Sigma index is %4.4f\n',sigma');
fprintf('The Var(G) is %4.4f\n',variance');
%fprintf('The IRM1(G) Index is %4.4f\n',IRM1');
% fprintf('The IR1(G) Index is %4.4f\n',IR1');
% fprintf('The IRM2(G) Index is %4.4f\n',IRM2');
% fprintf('The IR2(G) Index is %4.4f\n',IR2');
% fprintf('The IRC(G) Index is %4.4f\n',IRC');
% fprintf('The IRF(G) Index is %4.4f\n',IRF');
% fprintf('The IRFW(G) Index is %4.4f\n',IRFW');
% fprintf('The IRV(G) Index is %4.4f\n',IRV');
fprintf('The IRA(G) Index is %4.4f\n',IRA');
% fprintf('The IRB(G) Index is %4.4f\n',IRB');
% fprintf('The (G) Index is %4.4f\n',IRL');
fprintf('The IRDIF(G) Index is %4.4f\n',IRDIF');
fprintf('The IRZ(G) Index is %4.4f\n',IRZ');
fprintf('The IRLU(G) Index is %4.4f\n',IRLU');
fprintf('The IRLF(G) Index is %4.4f\n',IRLF');
fprintf('The IRLA(G) Index is %4.4f\n',IRLA');
fprintf('The IRD1(G) Index is %4.4f\n',IRD1');
fprintf('The IRGA(G) Index is %4.4f\n',IRGA');
fprintf('The IRT1(G) Index is %4.4f\n',IRT1');
fprintf('The IRT2(G) Index is %4.4f\n',IRT2');