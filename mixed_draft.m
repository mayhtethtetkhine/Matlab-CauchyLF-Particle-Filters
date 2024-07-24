%If our proposed better, less variance in ksdensity (narrower plot).
clear all
min_rl0 = 10;
gmin_rl0 = 10;
min_rl02 = 10;
gmin_rl02 = 10;
for bigtemploop=2:2 %vk variance 1 2 4, Line 13
for temploop=1:1 %no. of parts 100 200 300..., Line 7
    clear x1 w1 res res_new x2 w2 w2old wg2old
   % clear x1 xg1 w1 wg1 res res_new x2 w2 w2old
rand('seed',1); randn('seed',1); 
rloop0=500; %Use with bestidx rloop from rloop
stepnum=100;


%stepnum=200;

npart1=20*temploop; npart=20*temploop;
rnum(1:rloop0)=round(13000*rand(rloop0,1)); %seed number for each run
%xk = 0.5(xk-1)+(25(xk-1)/(1+(xk-1)^2))+8cos(1.2(k-1))+uk
%yk = 0.05(xk)^2+vk
%sigma_u=sqrt(2); %uk ~ N(0,2) and for x0i ~ N(0.1,2)
sigma_u=sqrt(0.1);
sigma_v=sqrt(0.1^(bigtemploop-1)); %vk ~ N(0,1).

%gamma = 0.8*(sigma_v);
gamma = sqrt(2/pi)* sigma_v;
%Clean states and measurements.
x0=0; %Prior

cn=[1 12 7]; %1 12 7
xc(1)=cn(1)*x0+((cn(2)*x0)/(1+(x0^2)))+cn(3);
xc
for cc=2:stepnum
xc(cc)=cn(1)*xc(cc-1)+((cn(2)*xc(cc-1))/(1+(xc(cc-1)^2)))+cn(3)*cos(1.2*(cc-1));
end
xc=xc+[zeros(40,1)' 30*ones(60,1)'];
%xc=xc+[zeros(100,1)' 30*ones(100,1)'];

%If all w=0, reset them to 1/npart and continue the
%simulation for the sake of comparison.
yc(1:stepnum)=0.05*xc(1:stepnum).^2;

%xn=awgn(xc,0,'measured'); %without measured, snr lower than the assigned
%noise=xn-xc;
%snr=10*log10((norm(xc)).^2/(norm(noise)).^2);
%xc=xn;

rloop=50; %Real loop numbers from initial loops for rand and randn

%copynp=zeros(npart,stepnum,rloop); %Number of replicas in resampling

%Start simulations
%tic %Start time here!
for rl0=1:rloop;
rand('seed',rnum(rl0)); randn('seed',rnum(rl0));

for cc=1:stepnum
    noise=randn(1,1)*sigma_v;
    yn(cc)=yc(cc)+noise;
    snr(cc)=20*log10(norm(yc(cc))/norm(noise)); %Clean signal to noise
end
avgsnr=mean(snr); %Fixed SNR but unknown

x(1:npart1,1)=x0+sigma_u*randn(npart1,1); %Initial parts x0i ~ N(0.1,2)
%std(1:npart1,1)=randn(npart1,1); %Randomly estimate variance.

xg(1:npart1,1)=x0+sigma_u*randn(npart1,1);
%Update parts/states
for i=1; [temploop sigma_v rl0 i] %this i means step 1
for ip=1:npart1
%Perturb x0i to get x1i.
x(ip,1)=cn(1)*x(ip,1)+((cn(2)*x(ip,1))/(1+(x(ip,1))^2))+cn(3)+sigma_u*randn(1,1);
xg(ip,1)=cn(1)*xg(ip,1)+((cn(2)*xg(ip,1))/(1+(xg(ip,1))^2))+cn(3)+sigma_u*randn(1,1);

%x1i=0.5*x0i+...8cos(1.2*0)
y(ip,1)=0.05*x(ip,1)^2; % y obtained from perturbed x
yg(ip,1)=0.05*xg(ip,1)^2;
%xxs=randn(1,1).^2; %N(0,1) %xxs=randn(length(yn(i)),1).^2;
%ssx=sum(xxs);
%ssx=ssx/sum((yn(i)-y(ip,1)).^2);
%std(ip,1)=1/sqrt(ssx);
%loglik(ip,1)=(-0.5/(std(ip,1)^2))*sum((yn(i)-y(ip,1)).^2);


%gaussian
logglik(ip,1)=(-0.5/(sigma_v^2))*sum((yn(i)-yg(ip,1)).^2);

%cauchy
loglik(ip,1) = -log(1+(sum(yn(i)-y(ip,1))/gamma).^2);

%loglik(ip,1) = log((1+((yn(i)-y(ip,1))/gamma).^2))^(-1);

%plik(ip,1) = loglik(ip,1);

end
[a,b]=max(loglik);
%gaussian
[ag, bg] = max(logglik);

for aa=1:npart1
   plik=exp(loglik-a); %plik=exp(loglik);
   pglik=exp(logglik-ag);
   
end

%cauchy
w=plik;

%gaussian
wg = pglik;

%------


res=[x w]; mapest(i,rl0)=res(b,1); %MAP just max w (No need normalization)

gres=[xg wg]; gmapest(i,rl0)=gres(bg,1);

%[res1,idx0]=sortrows(res,2); res1=flipud(res1); idx0=flipud(idx0); %Desc sort
%loglik=loglik(idx0(1:npart)); res1=res1(1:npart,:); res=res1; % Keep first npart parts for next states
%Non-normalized of first npart, used in A-R GA

%cauchy
w=res(:,2);
%gaussian
wg = gres(:,2);

%w=res(:,2)/npart; %w0=1/npart
res(:,2)=w;
oldsum=sum(w);

gres(:,2) = wg;
oldsumg=sum(wg);

%cauchy
%Neff and MMSE %Use norm w
%neff1(i,rl0)=(oldsum)^2/sum(w.^2); %Neff, Eq to 1/sum(normw.^2)
%wmest(i,rl0)=sum(res(:,1).*res(:,2)/oldsum); %sum(x(i)w(i))

% clear a aa b x y std plik; %clear old npart1 plik

%gaussian
%gneff1(i,rl0) = (oldsumg)^2/sum(wg.^2);
%gwmest(i,rl0)= sum(gres(:,1).*gres(:,2)/oldsumg);

clear a ag aa b bg x y yg std plik pglik; %clear old npart1 plik

x1(:,i,rl0)=res(:,1);
%w1(:,i,rl0)=res(:,2);
w1(:,i,rl0)=res(:,2)/oldsum;
varx1(i,rl0)=var(res(:,1)); 
%varw1(i,rl0)=var(res(:,2)); %Var of non-normalized w
varw1(i,rl0)=var(res(:,2)/oldsum);


xg1(:,i,rl0)=gres(:,1);
wg1(:,i,rl0)=gres(:,2)/oldsumg;
gvarx1(i,rl0)=var(gres(:,1));
gvarw1(i,rl0)=var(gres(:,2)/oldsumg);


%fprintf('%.6f',res1(tsnum,2))

%APF Ahwiadi 2020 %Use non-norm w in this process
[xbest,ii]=max(res(:,2)); clear xbest;
xbest=res(ii,1);
wbest=res(ii,2);
stdpart=std(res(:,1));  


%gaussian
[xgbest,ii]=max(gres(:,2)); clear xgbest;
xgbest=gres(ii,1);
wgbest=gres(ii,2);
stdpartg=std(gres(:,1));  


for ip=1:npart
    
  %-- cauchy
   if res(ip,2)<oldsum/npart
       %Find borders
       if res(ip,1)>=xbest
          lowb=xbest-stdpart;
          upb=res(ip,1)+stdpart;       
       else
          lowb=res(ip,1)-stdpart;
          upb=xbest+stdpart;       
       end
       temppl=res(ip,2);
       while temppl<oldsum/npart
          tempx=lowb+(upb-lowb)*rand(1,1);
          tempy=0.05*tempx^2;
          %templl=(-0.5/(sigma_v^2))*sum((yn(i)-tempy).^2);
            templl = -log(1+(sum(yn(i)-tempy)/gamma).^2);
          temppl=exp(templl-max(loglik)); %max loglik of prior APF for consistency.          
       end
       res(ip,:)=[tempx temppl];
       %Re-evaluate best particle and std of all particle values
       if res(ip,2)>wbest
           xbest=res(ip,1);
           wbest=res(ip,2);
       end
       stdpart=std(res(:,1));
   end
%-- 


%-- gaussian

 if gres(ip,2)<oldsumg/npart
       %Find borders
       if gres(ip,1)>=xgbest
          lowgb=xgbest-stdpartg;
          upgb=gres(ip,1)+stdpartg;       
       else
          lowgb=gres(ip,1)-stdpartg;
          upgb=xgbest+stdpartg;       
       end
       temppl=gres(ip,2);
       while temppl<oldsumg/npart
          tempx=lowgb+(upgb-lowgb)*rand(1,1);
          tempy=0.05*tempx^2;
          templl=(-0.5/(sigma_v^2))*sum((yn(i)-tempy).^2);
          temppl=exp(templl-max(logglik)); %max loglik of prior APF for consistency.          
       end
       gres(ip,:)=[tempx temppl];
       %Re-evaluate best particle and std of all particle values
       if gres(ip,2)>wgbest
           xgbest=gres(ip,1);
           wgbest=gres(ip,2);
       end
       stdpartg=std(gres(:,1));
   end


%--
    
end

    res_new=res; 
    gres_new = gres;


    clear idx0 xxs
    
%clear a idx1 c aa d res1 resh plik loglik y;
clear a ag idx1 c aa d res1 gres1 resh plik pglick loglik logglik y yg;

%------------------

x2(:,i,rl0)=res_new(:,1);
w2(:,i,rl0)=res_new(:,2); [a,b]=max(res_new(:,2));

w2old(:,i,rl0)=w2(:,i,rl0)/oldsum;

newsum=sum(res_new(:,2)); %Sum of all new weights

wmest2(i,rl0)=sum(res_new(:,1).*res_new(:,2)/newsum);
mapest2(i,rl0)=res_new(b,1); clear a b;
varx2(i,rl0)=var(res_new(:,1));
%neff2(i,rl0)=(newsum)^2/sum(res_new(:,2).^2); %Neff, Eq to 1/sum(normw.^2)

xmean(i,rl0)=sum(res_new(:,1))/npart;

%---- gaussian

xg2(:,i,rl0)=gres_new(:,1);
wg2(:,i,rl0)=gres_new(:,2); [ag,bg]=max(gres_new(:,2));

wg2old(:,i,rl0)=wg2(:,i,rl0)/oldsumg;

newsumg=sum(gres_new(:,2)); %Sum of all new weights

gwmest2(i,rl0)=sum(gres_new(:,1).*gres_new(:,2)/newsumg);
gmapest2(i,rl0)=gres_new(bg,1); clear ag bg;
gvarx2(i,rl0)=var(gres_new(:,1));
%gneff2(i,rl0)=(newsumg)^2/sum(gres_new(:,2).^2); %Neff, Eq to 1/sum(normw.^2)

gxmean(i,rl0)=sum(gres_new(:,1))/npart;


%---
 
%  [A, B] = ksdensity(res_new(:,1)); %Replot by plot(B,A) %Smoother than histogram
%  Aplot(:,rl0,i)=A; Bplot(:,rl0,i)=B;
%  %Eval the density est at 100 points covering data range.
%  % A: Vector of density values. % B: Set of 100 points.
%  [A1 B1]=max(A); modd2(i,rl0)=B(B1);

%Diversify particles before predicting next values
%-------------
%-- cauchy
stdpart2=std(res_new(:,1));

for ip=1:npart
    %res_new(ip,1)=2*res_new(ip,1)-xbest; %Opposite sign obtained    
    %res_new(ip,1)=res_new(ip,1)+randn(1,1)*stdpart2;
    res_new(ip,1)=res_new(ip,1)+randn(1,1)/stdpart2; %Use this line acc to paper.    
end

%-- gaussian
stdpart2g=std(gres_new(:,1));
for ip=1:npart
    %res_new(ip,1)=2*res_new(ip,1)-xbest; %Opposite sign obtained    
    %res_new(ip,1)=res_new(ip,1)+randn(1,1)*stdpart2;
    gres_new(ip,1)=gres_new(ip,1)+randn(1,1)/stdpart2g; %Use this line acc to paper.    
end

end

for i=2:stepnum; [temploop sigma_v rl0 i]
%clear plik loglik
clear plik pglik loglik logglik
for ip=1:npart

x(ip,1)=cn(1)*res_new(ip,1)+((cn(2)*res_new(ip,1))/(1+(res_new(ip,1))^2))+cn(3)*cos(1.2*(i-1))+sigma_u*randn(1,1);

%At i=2, x2i=0.5x1i+(25x1i/(1+x1i^2))+8cos(1.2*1)+noise
y(ip,1)=0.05*x(ip,1)^2; %y from perturbed particle x
%xxs=randn(1,1).^2; %N(0,1) %xxs=randn(length(yn(i)),1).^2;
%ssx=sum(xxs);
%ssx=ssx/sum((yn(i)-y(ip,1)).^2);
%std(ip,1)=1/sqrt(ssx);
%loglik(ip,1)=(-0.5/(std(ip,1)^2))*sum((yn(i)-y(ip,1)).^2);
%loglik(ip,1)=(-0.5/(sigma_v^2))*sum((yn(i)-y(ip,1)).^2);

loglik(ip,1) = -log(1+(sum(yn(i)-y(ip,1))/gamma).^2);


%--- gaussian

xg(ip,1)=cn(1)*gres_new(ip,1)+((cn(2)*gres_new(ip,1))/(1+(gres_new(ip,1))^2))+cn(3)*cos(1.2*(i-1))+sigma_u*randn(1,1);
yg(ip,1)=0.05*xg(ip,1)^2; %y from perturbed particle x
logglik(ip,1)=(-0.5/(sigma_v^2))*sum((yn(i)-yg(ip,1)).^2);

%---




end
%
[a,b]=max(loglik);
[ag,bg]=max(logglik);
for aa=1:npart
    plik=exp(loglik-a); 
   % plik=exp(loglik);

   %--gaussian

    pglik=exp(logglik-ag);
   %--
 
end



w=plik(:,1); %Reset wk-1 to 1/npart???
oldsum=sum(w); 
res(:,2)=w;
res=[x w]; %res_plot(:,:,i)=res;

%-- gaussian

wg = pglik(:,1);
oldsumg = sum(wg);
gres(:,2) = wg;
gres = [xg wg];


%--
%Neff, MAP, and MMSE
%neff1(i,rl0)=(oldsum)^2/sum(w.^2); %Eq to 1/sum(normw.^2)
mapest(i,rl0)=res(b,1); %Use index b for finding argmax
%wmest(i,rl0)=sum(res(:,1).*res(:,2)/oldsum); %sum(x(i)w(i))


%-- gaussian

%gneff1(i,rl0)=(oldsumg)^2/sum(wg.^2); %Eq to 1/sum(normw.^2)
gmapest(i,rl0)=gres(bg,1); %Use index b for finding argmax
%gwmest(i,rl0)=sum(gres(:,1).*gres(:,2)/oldsumg); %sum(x(i)w(i))

%----

%clear a aa b x std y;
clear a ag aa b bg xg std y;


x1(:,i,rl0)=res(:,1);
%w1(:,i,rl0)=res(:,2);
w1(:,i,rl0)=res(:,2)/oldsum;
varx1(i,rl0)=var(res(:,1));
%varw1(i,rl0)=var(res(:,2)); %Var of non-normalized w
varw1(i,rl0)=var(res(:,2)/oldsum);


%-- gaussian

xg1(:,i,rl0)=gres(:,1);
wg1(:,i,rl0)=gres(:,2)/oldsumg;
gvarx1(i,rl0)=var(gres(:,1));
gvarw1(i,rl0)=var(gres(:,2)/oldsumg);

%----

%APF Ahwiadi 2020 %Use non-norm w in this process
[xbest,ii]=max(res(:,2)); clear xbest;
xbest=res(ii,1);
wbest=res(ii,2);
stdpart=std(res(:,1));  

%-- gaussian

[xgbest,ii]=max(gres(:,2)); clear xgbest;
xgbest=gres(ii,1);
wgbest=gres(ii,2);
stdpartg=std(gres(:,1));  

%------


for ip=1:npart
  %---cauchy
   if res(ip,2)<oldsum/npart
       %Find borders
       if res(ip,1)>=xbest
          lowb=xbest-stdpart;
          upb=res(ip,1)+stdpart;       
       else
          lowb=res(ip,1)-stdpart;
          upb=xbest+stdpart;       
       end
       temppl=res(ip,2);
       while temppl<oldsum/npart
          tempx=lowb+(upb-lowb)*rand(1,1);
          tempy=0.05*tempx^2;
         % templl=(-0.5/(sigma_v^2))*sum((yn(i)-tempy).^2);
           templl = -log(1+(sum(yn(i)-tempy)/gamma).^2);
          temppl=exp(templl-max(loglik)); %max loglik of prior APF for consistency.          
       end
       res(ip,:)=[tempx temppl];
       %Re-evaluate best particle and std of all particle values
       if res(ip,2)>wbest
           xbest=res(ip,1);
           wbest=res(ip,2);
       end
       stdpart=std(res(:,1));
   end
%--- 

%-- gaussian

 if gres(ip,2)<oldsumg/npart
       %Find borders
       if gres(ip,1)>=xgbest
          lowgb=xgbest-stdpartg;
          upgb=gres(ip,1)+stdpartg;       
       else
          lowgb=gres(ip,1)-stdpartg;
          upgb=xgbest+stdpartg;       
       end
       temppl=gres(ip,2);
       while temppl<oldsumg/npart
          tempx=lowgb+(upgb-lowgb)*rand(1,1);
          tempy=0.05*tempx^2;
          templl=(-0.5/(sigma_v^2))*sum((yn(i)-tempy).^2);
          % templl = -log(1+(sum(yn(i)-tempy)/gamma).^2);
          temppl=exp(templl-max(logglik)); %max loglik of prior APF for consistency.          
       end
       gres(ip,:)=[tempx temppl];
       %Re-evaluate best particle and std of all particle values
       if gres(ip,2)>wgbest
           xgbest=gres(ip,1);
           wgbest=gres(ip,2);
       end
       stdpartg=std(gres(:,1));
   end


%--


end
    res_new=res;   
   gres_new = gres;

    clear idx0 xxs
    
%clear a idx1 c aa d res1 resh plik loglik y;
clear a ag idx1 c aa d res1 resh plik pglik loglik logglik y yg;

x2(:,i,rl0)=res_new(:,1);
w2(:,i,rl0)=res_new(:,2); [a,b]=max(res_new(:,2));

w2old(:,i,rl0)=w2(:,i,rl0)/oldsum;

newsum=sum(res_new(:,2)); %Sum of all new weights
wmest2(i,rl0)=sum(res_new(:,1).*res_new(:,2)/newsum);
mapest2(i,rl0)=res_new(b,1); clear a b;
varx2(i,rl0)=var(res_new(:,1));
%neff2(i,rl0)=(newsum)^2/sum(res_new(:,2).^2); %Neff, Eq to 1/sum(normw.^2)

xmean(i,rl0)=sum(res_new(:,1))/npart;

%-- gaussian

xg2(:,i,rl0)=gres_new(:,1);
wg2(:,i,rl0)=gres_new(:,2); [ag,bg]=max(gres_new(:,2));
wg2old(:,i,rl0)=wg2(:,i,rl0)/oldsumg;
newsumg=sum(gres_new(:,2)); %Sum of all new weights

gwmest2(i,rl0)=sum(gres_new(:,1).*gres_new(:,2)/newsumg);
gmapest2(i,rl0)=gres_new(bg,1); clear ag bg;
gvarx2(i,rl0)=var(gres_new(:,1));
%gneff2(i,rl0)=(newsumg)^2/sum(gres_new(:,2).^2); %Neff, Eq to 1/sum(normw.^2)

gxmean(i,rl0)=sum(gres_new(:,1))/npart;

%---

 
%  [A, B] = ksdensity(res_new(:,1)); %Replot by plot(B,A) %Smoother than histogram
%  Aplot(:,rl0,i)=A; Bplot(:,rl0,i)=B;
%  %Eval the density est at 100 points covering data range.
%  % A: Vector of density values. % B: Set of 100 points.
%  [A1 B1]=max(A); modd2(i,rl0)=B(B1);

%Diversify particles before predicting next values
%--cauchy
stdpart2=std(res_new(:,1));
for ip=1:npart
    res_new(ip,1)=res_new(ip,1)+randn(1,1)/stdpart2; %Use this line acc to paper.    
end

%-- gaussian
stdpart2g=std(gres_new(:,1));
for ip=1:npart
    gres_new(ip,1)=gres_new(ip,1)+randn(1,1)/stdpart2g; %Use this line acc to paper.    
end

end

rmse1(rl0)=sqrt(immse(mapest2(:,rl0),xc.'));

if rl0 ~= 1
if rmse1(rl0)<min_rl0(1,1)
        min_rl0 = [rmse1(rl0) rl0];
end
end


abserr1(rl0)=sum(abs(mapest2(:,rl0)-xc.'))/stepnum;

if rl0 ~= 1
if abserr1(rl0)<min_rl02(1,1)
        min_rl02 = [abserr1(rl0) rl0];
end
end


%--gaussian
grmse1(rl0)=sqrt(immse(gmapest2(:,rl0),xc.'));

if rl0 ~= 1
if grmse1(rl0)<gmin_rl0(1,1)
        gmin_rl0 = [grmse1(rl0) rl0];
end
end

gabserr1(rl0)=sum(abs(gmapest2(:,rl0)-xc.'))/stepnum;

if rl0 ~= 1
if gabserr1(rl0)<gmin_rl02(1,1)
        gmin_rl02 = [gabserr1(rl0) rl0];
end
end

%--
end
%toc


%-----need --------


figure('Name', 'Cauchy');%0-100

for rl0=1:rloop
    if rl0 == 1 
    plot(0:stepnum,1*([x0 (mapest2(:,rl0)).']).','r', 'DisplayName', 'Estimations', 'Parent', gca);
    
    else 
    plot(0:stepnum,1*([x0 (mapest2(:,rl0)).']).','r');
    
    end

hold on;
end

plot(0:stepnum,([x0 xc]).','k-o','LineWidth',1, 'DisplayName', 'True Value' ,'Parent', gca); hold on;

title('MAP using Cauchy Likelihood Function');

legendEntries = findobj('Type', 'line', 'Parent', gca, '-and', '-not', 'DisplayName', ''); % Find lines without DisplayName
legend(legendEntries, 'Location', 'Best'); 
%legend('show', 'AutoUpdate', 'off');
fontsize(gca,30,"pixels")
xlim([0 100]);
ylim([-30 100]);
xlabel('Time-step (k)');
ylabel('x_k');

%-- gaussian
figure('Name', 'Gaussian'); %0-100
for rl0=1:rloop
   if rl0 == 1

plot(0:stepnum,1*([x0 (gmapest2(:,rl0)).']).','r', 'DisplayName', 'Estimations');
   else
plot(0:stepnum,1*([x0 (gmapest2(:,rl0)).']).','r');
   end
hold on;
end

plot(0:stepnum,([x0 xc]).','k-o','LineWidth',1,'DisplayName', 'True Value'); hold on;
title('MAP using Gaussian Likelihood Function');

secondFigureAxes = gca;

%legendEntries2 = legendEntries; 
legendEntries2 =  findobj('Type', 'line', 'Parent', secondFigureAxes, '-and', '-not', 'DisplayName', ''); 
legend(legendEntries2, 'Location', 'best'); 

fontsize(gca,30,"pixels")
xlim([0 100]);
ylim([-30 100]);
xlabel('Time-step (k)');
ylabel('x_k');
%--




figure('Name', 'Minimum RMSE comparison'); %0-100

plot(0:stepnum,([x0 xc]).','k','LineWidth',1.5); 
hold on;
plot(0:stepnum,1*([x0 (mapest2(:,min_rl0(1,2))).']).','g--','LineWidth',1);
hold on;
plot(0:stepnum,1*([x0 (gmapest2(:,gmin_rl0(1,2))).']).','r:','LineWidth',1);
hold on;

title('Minimum RMSE', FontSize= 20);
legend('True', 'CLF', 'GLF');
fontsize(gca,30,"pixels")
xlabel('Time-step (k)', FontSize=30);
ylabel('x_k', FontSize=30);

%--

figure('Name', 'Minimum AbsErr comparison'); %0-100

plot(0:stepnum,([x0 xc]).','k','LineWidth',1); 
hold on;
plot(0:stepnum,1*([x0 (mapest2(:,min_rl02(1,2))).']).','g--','LineWidth',1);
hold on;
plot(0:stepnum,1*([x0 (gmapest2(:,gmin_rl02(1,2))).']).','r:','LineWidth',1);
hold on;

title('Minimum Absolute Error', FontSize= 20);
legend('True', 'CLF', 'GLF');
fontsize(gca,30,"pixels")
xlabel('Time-step (k)', FontSize= 30);
ylabel('x_k', FontSize= 30);


%-- needed

%}

rmse1

minrmse1=min(rmse1);     avgrmse1=mean(rmse1);     maxrmse1=max(rmse1); varrmse1=var(rmse1);
gminrmse1=min(grmse1);     gavgrmse1=mean(grmse1);     gmaxrmse1=max(grmse1); gvarrmse1=var(grmse1);
minabserr1=min(abserr1); avgabserr1=mean(abserr1); maxabserr1=max(abserr1); varabserr1=var(abserr1);
gminabserr1=min(gabserr1); gavgabserr1=mean(gabserr1); gmaxabserr1=max(gabserr1); gvarabserr1=var(gabserr1);

end
end

%{
[ min_rmse_of_c  avg_rmse_of_c  max_rmse_of_c  variance_rmse_of_c
min_rmse_of_gau  avg_rmse_of_gau  max_rmse_of_gau  variance_rmse_of_gau
min_abserr_of_c  avg_abserr_of_c  max_abserr_of_c  variance_abserr_of_c
min_abser_of_gau  avg_abser_of_gau  max_abser_of_gau variance_abser_of_gau]

%}

[minrmse1 avgrmse1 maxrmse1 varrmse1;...
gminrmse1 gavgrmse1 gmaxrmse1 gvarrmse1;...
minabserr1 avgabserr1 maxabserr1 varabserr1;...
gminabserr1 gavgabserr1 gmaxabserr1 gvarabserr1]

%avg of CLF - avg of GLF
%[avgrmse1-gavgrmse1  avgabserr1-gavgabserr1]

[avgrmse1 gavgrmse1 avgabserr1 gavgabserr1; ...
    minrmse1 gminrmse1 minabserr1 gminabserr1]

%{

sigma = sqrt(0.2)
CLF avg rmse and abserr < GLF ###
but variance of CLF > GLF

sigma = sqrt(3)
CLF avg rmse and abserr > GLF
but variance of CLF > GLF
  
%}
