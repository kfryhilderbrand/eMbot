% Define a set of 64-element vectors that represent all possisle
% connections of the two valve settings on six cylinders:
% [0, A]: 
%number of cylinders
k=6;
%number of states
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

% Define the different cross sections for the six cylinders

%% IROS values 
A1 = 361;
A2 = -676;
A3 = 991;
A4 = 361;
A5 = -676;
A6 = 991;

%% Close to J=0
% A1 = 361;
% A2 = 676;
% A3 = -991;
% A4 = 361;
% A5 = 676;
% A6 = -991;

%%
[~, n]=size(s1);
T=zeros(1,n);
%T stores all possible transmission values
for i=1:n
    T(i)= -(A4*s4(i) + A5*s5(i) + A6*s6(i))/(A1*s1(i) + A2*s2(i) + A3*s3(i));
end

%Alpha_Beta stores all finite and unique transmission values
Alpha_Beta=unique(T);
Alpha_Beta=sort(Alpha_Beta(isfinite(Alpha_Beta)));
disp(vpa(Alpha_Beta,3))

%% Number of Unique Transmissions
disp('Number of independent transmissions:')
[m, n] = size(Alpha_Beta);
disp(n)
Beta_Alpha=sort(1./Alpha_Beta);

%% Minimum and Maximum Transmission
disp('Maximum Transmission:')
disp(vpa(max(Alpha_Beta),3))
disp('Minimum Transmission:')
disp(vpa(min(Alpha_Beta),3))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% clf
% hold on
% grid on
% bar(T,'r')
% xlabel('Combination #')
% ylabel('Transmission ratio T')
% 
% figure(2)
% clf
% hold on
% grid on
% bar(Alpha_Beta,'b')
% xlabel('Sorted Valve Setting Combination')
% ylabel('Transmission ratio T')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% figure(1)
% clf
% hold on
% grid on
% bar(ang*180/pi,'b')
% axis([0,50,-200,200])
% xlabel('Combination #')
% ylabel('joint space angle')
% 
% figure(2)
% clf
% hold on
% grid on
% bar(sort(tan(ang)),'r')
% %bar(tan(ang),'r')
% axis([0,50,-30,30])
% xlabel('Combination #')
% ylabel('Transmission ratio T')
% 
% figure(3)
% clf
% hold on
% grid on
% axis equal
% for i = 1:length(ang)
%     plot([-cos(ang(i)),cos(ang(i))],[-sin(ang(i)),sin(ang(i))]);
% end
% xlabel('\alpha_{in}')
% ylabel('\alpha_{out}')
% axis([-2,2,-2,2]);
