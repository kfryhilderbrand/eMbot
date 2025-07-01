% Define a set of 64-element vectors that represent all possisle
% connections of the two valve settings on six cylinders:
% [0, A]: 
s = [0, 1];
s1 = repmat(s,2^0,2^5);
s1 = reshape(s1, 2^0, 2^6);
s2 = repmat(s,2^1,2^4);
s2 = reshape(s2, 2^0, 2^6);
s3 = repmat(s,2^2,2^3);
s3 = reshape(s3, 2^0, 2^6);
s4 = repmat(s,2^3,2^2);
s4 = reshape(s4, 2^0, 2^6);
s5 = repmat(s,2^4,2^1);
s5 = reshape(s5, 2^0, 2^6);
s6 = repmat(s,2^5,2^0);
s6 = reshape(s6, 2^0, 2^6);


%% %% General Case
% % Define the different cross sections for the six cylinders
syms j1 j2 j3 j4 j5 j6;
[~, n]=size(s1);
T=zeros(1,n);
T=sym(T);
for i=1:n
T(i)= -(j4*s4(i) + j5*s5(i) + j6*s6(i))/(j1*s1(i) + j2*s2(i) + j3*s3(i));
    if((j4*s4(i) + j5*s5(i) + j6*s6(i)==0)&&((j1*s1(i) + j2*s2(i) + j3*s3(i)==0)))
        T(i)=Inf;
    end
end
sort(T);
Alpha_Beta=unique(T);
disp(Alpha_Beta)
disp('Number of independent transmissions (including INF):')
[~, n] = size(Alpha_Beta);
disp(n)
Beta_Alpha=sort(1./Alpha_Beta);

%% Lever arms sum to zero
% Define the different cross sections for the six cylinders
syms j1 j2 j3 j4 j5 j6;
j3 = -(j1+j2);
j6 = -(j4+j5);
[~, n]=size(s1);
T=zeros(1,n);
T=sym(T);
for i=1:n
T(i)= -(j4*s4(i) + j5*s5(i) + j6*s6(i))/(j1*s1(i) + j2*s2(i) + j3*s3(i));
    if((j4*s4(i) + j5*s5(i) + j6*s6(i)==0)&&((j1*s1(i) + j2*s2(i) + j3*s3(i)==0)))
        T(i)=Inf;
    end
end
sort(T);
Alpha_Beta=unique(T);
disp(Alpha_Beta)
disp('Number of independent transmissions (including INF):')
[~, n] = size(Alpha_Beta);
disp(n)
Beta_Alpha=sort(1./Alpha_Beta);


%% Identical Lever arms sum to zero
% Define the different cross sections for the six cylinders
syms j1 j2 j3 j4 j5 j6;
j3 = -(j1+j2);
j4 = j1;
j5 = j2;
j6 = -(j1+j2);
[~, n]=size(s1);
T=zeros(1,n);
T=sym(T);
for i=1:n
T(i)= -(j4*s4(i) + j5*s5(i) + j6*s6(i))/(j1*s1(i) + j2*s2(i) + j3*s3(i));
    if((j4*s4(i) + j5*s5(i) + j6*s6(i)==0)&&((j1*s1(i) + j2*s2(i) + j3*s3(i)==0)))
        T(i)=Inf;
    end
end
sort(T);
Alpha_Beta=unique(T);
disp(Alpha_Beta)
disp('Number of independent transmissions (including INF):')
[~, n] = size(Alpha_Beta);
disp(n)
Beta_Alpha=sort(1./Alpha_Beta);

%% Fully Constrained
% % Define the different cross sections for the six cylinders
syms j1 j2 j3 j4 j5 j6;
j1 = 1;
j3 = -(j1+j2);
j4 = 1;
j5 = j2;
j6 = -(j1+j2);
[~, n]=size(s1);
T=zeros(1,n);
T=sym(T);
for i=1:n
T(i)= -(j4*s4(i) + j5*s5(i) + j6*s6(i))/(j1*s1(i) + j2*s2(i) + j3*s3(i));
    if((j4*s4(i) + j5*s5(i) + j6*s6(i)==0)&&((j1*s1(i) + j2*s2(i) + j3*s3(i)==0)))
        T(i)=Inf;
    end
end
sort(T);
Alpha_Beta=unique(T);
disp(Alpha_Beta)
disp('Number of independent transmissions (including INF):')
[~, n] = size(Alpha_Beta);
disp(n)
Beta_Alpha=sort(1./Alpha_Beta);


%%
%Positive Plots
figure(1)
j2=0.1:0.1:5;
plot(j2,subs(Alpha_Beta(3)),'b',j2,subs(Alpha_Beta(4)),'k',j2,subs(Alpha_Beta(8)),'y',...
j2,subs(Alpha_Beta(11)),'c',j2,subs(Alpha_Beta(12)),'r',j2,subs(Alpha_Beta(13)),'g',...
j2,subs(Alpha_Beta(14)),'m')
legend('1','j2','1/(j2+1)','j2+1','1/j2','j2/(j2+1)','(j2+1)/j2')
title('Positive Transmission Trajectories')
xlabel('Choice of j2')
ylabel('Possible Transmission Ratios')


%Negative Plots
figure(2)
j2=0.1:0.1:5;
plot(j2,subs(Alpha_Beta(1)),'-b',j2,subs(Alpha_Beta(6)),'-k',j2,subs(Alpha_Beta(7)),'-y',...
j2,subs(Alpha_Beta(9)),'-c',j2,subs(Alpha_Beta(10)),'-r',j2,subs(Alpha_Beta(15)),'-g',...
j2,subs(Alpha_Beta(16)),'-m')
hleg=legend('-1','-j2-1','-1/j2','-1/(j2+1)','-j2','-j2/(j2+1)','-(j2+1)/j2');
set(hleg,'Location','SouthEast')
title('Negative Transmission Trajectories')
xlabel('Choice of j2')
ylabel('Possible Transmission Ratios')

