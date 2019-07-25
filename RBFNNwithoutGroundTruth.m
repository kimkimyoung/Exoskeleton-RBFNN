clear all;
close all;

encoder = importdata('encoder.txt');
imu_data = importdata('imu.txt');
imu_data = imu_data.data;
strain_gauge_data = importdata('strain_gauge.txt');

strain_gauge_data(:,2:4) = [];
encoder = encoder.';
[NUM,~] = size(encoder);
encoder_data = [];
sum = 0;
j = 1;
for i = 1:1:NUM
    
    encoder(i) = encoder(i) * 360 / 4000;
    sum = sum + encoder(i);
    if ( mod(i,10) == 0 )
        encoder_data(j) = sum / 10;
        j = j+1;
        sum = 0;
    end

end
encoder_data = encoder_data.';
encoder_data(1:133,:) = [];
encoder_data(641:end,:) = [];
imu_data(:,1:6) = [];
imu_data(:,2:7) = [];

strain_gauge_data(1:418,:) = [];
imu_data(1:460,:) = [];
% 
imu_data(641:end,:) = [];
strain_gauge_data(641:end,:) = [];
strain_gauge_data = (strain_gauge_data - 477)/5  ;

[row,~] = size(imu_data);
imu_data_vel = [];
alfa = 0.3;
for k = 1:1:row
    if k == 1
        imu_data_vel(k) = 0;
    else
        imu_data_vel(k) = alfa * imu_data(k) - ( 1-alfa ) * imu_data(k-1) ;
    end
end
imu_data_vel = imu_data_vel.';
        
% figure(1)
% s_plot = plot(strain_gauge_data);
% hold on;
% % figure(2)
% eplot = plot(encoder_data);
% hold on;
% iplot = plot(imu_data);
% hold on;
% legend([s_plot,eplot,iplot],'strain','encoder','imu');
% 
% figure(3)
% plot(imu_data_vel);

%%
% strain_gauge_data = strain_gauge_data / 400;
Dataset = [strain_gauge_data,imu_data,imu_data_vel];
Xd = encoder_data;
[row,col] = size(Dataset);

alfa = 0.1; % 动量因子
xite = 0.3; % 学习效率

% 初始化参数
n_center_vector = 10;
c = -2 + 4.* rand(col,n_center_vector);% [-2 2]
b = 3 * ones(n_center_vector,1);
w = 0.1 * ones(1, n_center_vector);
h = zeros(n_center_vector,1);

w_1 = w; w_2 = w_1;
c_1 = c; c_2 = c_1;
b_1 = b; b_2 = b_1;
d_w = 0 * w;
d_b = 0 * b;
d_c = 0 * c;

% loop without groud truth
for s = 1:1:row
    
    x = Dataset(s,:);
    if(s ~= 1)
        err = Dataset(s,2) - yout(s-1);
        eout = 0.5 * err^2;
        es(s) = eout;
        
        d_w = xite * h' * err;
        for j = 1:1:n_center_vector
            d_b(j) = xite * err * w(j) * h(j)' * (b(j)^-3) * norm(x.' - c(:,j))^2;
            for i = 1:1:col
                d_c(i,j) = xite * err * w(j) * h(j) * (x(i) - c(i,j)) * (b(j)^-2);
            end
        end
        w = w_1 + d_w + alfa * (w_1 - w_2);
        b = b_1 + d_b + alfa * (b_1 - b_2);
        c = c_1 + d_c + alfa * (c_1 - c_2);

        w_2 = w_1;
        w_1 = w;

        c_2 = c_1;
        c_1 = c;

        b_2 = b_1;
        b_1 = b;  
    end
    
    for j = 1:1:n_center_vector
        h(j) = exp( -norm( x.' - c(:,j))^2 / (2 * b(j)^2));
    end
    % output layer
    yout(s) = w * h;
end

figure(1)

m = plot(yout);
hold on;
i = plot(imu_data);
hold on;
legend([i,m],'imu','RBFNN');
xlabel('data number');
ylabel('angle(degree)');

figure(2)
plot(es);
xlabel('error');