clear all
sigma = [];


for i=0.01:0.01:0.2
    sigma(end+1)=sqrt(i);
end


avg_rmse_clf = [3.9358 2.9674 2.5282 2.2514 2.4748 2.4640 2.1250 2.5907 2.3218 2.5248 2.1872 2.2684 2.4764 2.6254 2.3270 2.1835 2.5995 2.6870 2.7770 2.5509];
avg_rmse_glf = [13.1468 6.2149 3.1898 2.8592 2.9670 2.2490 2.4408  2.4751 2.6281 2.4930 2.5489 2.5458 2.6946 2.4121 2.7485 2.7875 2.8245 2.7459 2.8089 2.9006];
avg_abserr_clf = [0.87 0.7009 0.5944 0.4970 0.6052 0.5482 0.4921 0.6111 0.5367 0.5961 0.5083 0.5371 0.688 0.6134 0.5515 0.5171 0.6271 0.6337 0.6814 0.6169];
avg_abserr_glf = [3.4878 2.7008 0.7060 0.6684 0.6916 0.5255 0.5620 0.5668 0.6183 0.5666 0.5826 0.5887 0.6173 0.5415 0.6303 0.6414 0.6452 0.6283 0.6572 0.6738];

figure('Name', 'Varying Sigma Values AVG Errors')
subplot(2,1,1);
plot(sigma,avg_rmse_clf,'k:*', 'LineWidth',1);
hold on;
plot(sigma,avg_rmse_glf,'r:o', 'LineWidth',1);
legend('RMSE of CLF', 'RMSE of GLF');
xlabel('Standard Deviation (Sigma)');
ylabel('Mean RMSE ');
title('(a) Average RMSE');

subplot(2,1,2);
plot(sigma,avg_abserr_clf,'k:*', 'LineWidth',1);
hold on;
plot(sigma,avg_abserr_glf,'r:o', 'LineWidth',1);
legend( 'AbsErr of CLF',  'AbsErr of GLF');
xlabel('Standard Deviation (Sigma)');
ylabel('Mean MAE ');
title('(b) Average MAE');


min_rmse_clf = [0.1225 0.1294 0.1104 0.1546 0.1087 0.1309 0.1487 0.16 0.15 0.2126 0.1576 0.17 0.1625 0.1715 0.21 0.2522 0.2473 0.3251 0.3499 0.2062];
min_rmse_glf = [0.38 0.2374 0.3503 0.2636 0.3659 0.2404 0.2502 0.3883 0.4357 0.7168 0.6224 0.8763 0.9535 0.2353 0.8551 1.0166 1.0524 1.1152 1.1773 1.4527];

min_abserr_clf = [0.067 0.0702 0.0713 0.0878 0.0666 0.0828 0.0872 0.0963 0.0935 0.1137 0.1024 0.1027 0.1072 0.1092 0.1301 0.1317 0.1447 0.1487 0.1425 0.1280];
min_abserr_glf = [0.1404 0.1011 0.1293 0.1123 0.1321 0.0907 0.1074 0.1545 0.1579 0.1719 0.1572 0.1823 0.2029 0.1142 0.2077 0.2026 0.2305 0.2225 0.2273 0.2896];

figure('Name', 'Varying Sigma Values MIN Errors');
subplot(2,1,1);
plot(sigma,min_rmse_clf,'k:*', 'LineWidth',1);
hold on;
plot(sigma,min_rmse_glf,'r:o', 'LineWidth',1);
legend('RMSE of CLF', 'RMSE of GLF');
xlabel('Standard Deviation (Sigma)');
ylabel('Minimum RMSE ');
title('(a) Minimum RMSE');


subplot(2,1,2);
plot(sigma,min_abserr_clf,'k:*', 'LineWidth',1);
hold on;
plot(sigma,min_abserr_glf,'r:o', 'LineWidth',1);
legend( 'AbsErr of CLF',  'AbsErr of GLF');
xlabel('Standard Deviation (Sigma)');
ylabel('Minimum MAE ');
title('(b) Minimum MAE');


particle = [];
for k=1000:1000:25000
particle(end+1) = k;
end



pavg_rmse_clf = [2.3294 2.3042 2.3190 2.1696 2.2234 2.6649 2.6147 2.7847 2.7642 2.5248 2.6902 2.6277 2.8373 2.7548 2.2318 2.4778 2.7298 2.4873 2.7221 2.2033 2.3656 2.2743 2.3648 2.5919 2.2785];
pavg_rmse_glf = [2.2938 2.3610 2.3287 2.4727 2.4646 2.5641 2.4490 2.6736 2.5034 2.4930 2.4271 2.4191 2.4908 2.5331 2.6097 1.9324 2.4479 2.4940 2.6028 2.3411 2.3339 2.5038 2.3896 2.2060 2.3862];

pavg_abserr_clf=[0.5320 0.5263 0.5108 0.4915 0.5091 0.6433 0.6407 0.6648 0.6651 0.5961 0.6533 0.6271 0.6866 0.6603 0.5257 0.5838 0.6373 0.5867 0.6378 0.5089 0.5586 0.5049 0.5521 0.6170 0.5296];
pavg_abserr_glf=[0.5175 0.5385 0.5322 0.5619 0.5719 0.5937 0.5553 0.6271 0.5710 0.5666 0.5596 0.5559 0.5696 0.5828 0.5945 0.4262 0.5584 0.5722 0.6088 0.5217 0.5298 0.5696 0.5377 0.4954 0.5401];

%{

figure('Name', 'Varying Particle Number')
subplot(2,2,1);
plot(particle,pavg_rmse_clf,'k:*', 'LineWidth',1);
hold on;
plot(particle,pavg_rmse_glf,'r:*', 'LineWidth',1);
legend('RMSE of CLF', 'RMSE of GLF');
xlabel('Number of Particle');
ylabel('Mean RMSE ');
title('Comparison of Average RMSE of Cauchy and Gaussian');

subplot(2,2,2);
plot(particle,pavg_abserr_clf,'k:*', 'LineWidth',1);
hold on;
plot(particle,pavg_abserr_glf,'r:*', 'LineWidth',1);
legend( 'AbsErr of CLF',  'AbsErr of GLF');
xlabel('Number of Particle');
ylabel('Mean AbsErr ');
title('Comparison of Average Absolute Error of Cauchy and Gaussian');

%}

pmin_rmse_clf = [0.2730 0.1457 0.2785 0.1574 0.1359 0.1526 0.1529 0.1877 0.1587 0.2126 0.1906 0.1730 0.2200 0.2262 0.1560 0.1701 0.1535 0.2260 0.1708 0.1576 0.1322 0.2275 0.1715 0.2382 0.2060];
pmin_rmse_glf = [0.4749 0.5356 0.5299 0.4614 0.4614 0.2145 0.7325 0.3267 0.7850 0.7168 0.2604 0.2307 0.7174 0.7463 0.8423 0.1471 0.2258 0.5057 0.4593 0.5355 0.2941 0.6867 0.3514 0.2963 0.2939];

pmin_abserr_clf = [0.1203 0.0944 0.1255 0.0941 0.0921 0.0991 0.0984 0.1069 0.1008 0.1137 0.1087 0.1056 0.1104 0.1284 0.0985 0.0960 0.0938 0.1241 0.1077 0.0965 0.0873 0.1054 0.0983 0.1124 0.1201];
pmin_abserr_glf = [0.1412 0.1624 0.1558 0.1806 0.1502 0.1172 0.2104 0.1294 0.2040 0.1719 0.1034 0.1159 0.1598 0.1673 0.1942 0.0875 0.1077 0.1538 0.1610 0.1575 0.1310 0.1793 0.1238 0.1256 0.1032]; 

figure('Name', 'Varying Particle Number MIN Errors');
subplot(2,1,1);
plot(particle,pmin_rmse_clf,'k:*', 'LineWidth',1);
hold on;
plot(particle,pmin_rmse_glf,'r:o', 'LineWidth',1);
legend('RMSE of CLF', 'RMSE of GLF');
xlabel('Number of Particle');
ylabel('Minimum RMSE ');
title('(a) Minimum RMSE');


subplot(2,1,2);
plot(particle,pmin_abserr_clf,'k:*', 'LineWidth',1);
hold on;
plot(particle,pmin_abserr_glf,'r:o', 'LineWidth',1);
legend( 'AbsErr of CLF',  'AbsErr of GLF');
xlabel('Number of Particle');
ylabel('Minimum MAE ');
title('(b) Minimum MAE');

%{
figure;
x = [-5:.1:5];
y = normpdf(x,0,1);
plot(x,y, 'r--')
hold on;

pd = makedist('tLocationScale','mu',0,'sigma',0.8, 'nu', 1);
x = -5:.1:5;
y = pdf(pd,x);
plot(x,y,'b','LineWidth',1)
legend ('Gaussian', 'Cauchy');
xlabel('x');
ylabel('f(x)');
%}

       