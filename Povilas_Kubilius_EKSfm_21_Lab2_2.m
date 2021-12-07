%Povilas_Kubilius_EKSfm-21
%Lab2_papildoma

clc
clear

x1 = linspace(0,1,25);
x2 = linspace(0,1,25);

x = [x1 x2];

d = ((1+ 0.6*sin(2/pi*x/0.7)) + 0.3*sin(2*pi*x))/2;
ERROR = 0;

%Atsitiktinis koeficientui generavimas

% Pirmasis sluoksnis
w11 = randn(1);
w12 = randn(1);

w21 = randn(1);
w22 = randn(1);

w31 = randn(1);
w32 = randn(1);

w41 = randn(1);
w42 = randn(1);

b1 = [randn(1) randn(1) randn(1) randn(1)];

%Antrasis sluoksnis
w11_2 = randn(1);
w12_2 = randn(1);
w13_2 = randn(1);
w14_2 = randn(1);

b1_2 = randn(1); 

%Masyvo dydzio nustatymas kodo efektyvumui pagerinti

v11 = zeros(1,length(x));
v12 = zeros(1,length(x));
v13 = zeros(1,length(x));
v14 = zeros(1,length(x));
v15 = zeros(1,length(x));

y11 = zeros(1,length(x));
y12 = zeros(1,length(x));
y13 = zeros(1,length(x));
y14 = zeros(1,length(x));
y15 = zeros(1,length(x));

v21 = zeros(1,length(x));

y   = zeros(1,length(x));

e   = zeros(1,length(x));


%Isejimo skaiciavimas
for i = 1:1:length(x)
    
         v11(i) = x(i) * w11 + x(i) * w12 + b1(1);
         v12(i) = x(i) * w21 + x(i) * w22 + b1(2); 
         v13(i) = x(i) * w31 + x(i) * w32 + b1(3); 
         v14(i) = x(i) * w41 + x(i) * w42 + b1(4); 
    
         %Pirmojo sluoksnio isejimai  
         y11(i) = 1/(1+exp(-v11(i)));
         y12(i) = 1/(1+exp(-v12(i)));
         y13(i) = 1/(1+exp(-v13(i)));
         y14(i) = 1/(1+exp(-v14(i)));

    v21(i) = (y11(i) * w11_2) + (y12(i) * w12_2) + (y13(i) * w13_2) + (y14(i) * w14_2) + b1_2;
    %Antrojo sluoksnio isejimas
    y(i) = v21(i); 
      
    e(i) = d(i) - y(i);
    ERROR = ERROR + abs(e(i));   
end

%mokymo zingsnis
h = 0.01; 

%Masyvo dydzio nustatymas kodo efektyvumui pagerinti

delta_out = zeros(1,length(x));

delta1 = zeros(1,length(x));
delta2 = zeros(1,length(x));
delta3 = zeros(1,length(x));
delta4 = zeros(1,length(x));
delta5 = zeros(1,length(x));

while ERROR > 0.15
    ERROR = 0;
    for i = 1:1:length(x)
        
        
        delta_out(i) = 1 * e(i);
        
        %Pasleptojo sluoksnio gradientas
        delta1(i) = y11(i) * ( 1 - y11(i) ) * ( delta_out(i) * w11_2 );
        delta2(i) = y12(i) * ( 1 - y12(i) ) * ( delta_out(i) * w12_2 );
        delta3(i) = y13(i) * ( 1 - y13(i) ) * ( delta_out(i) * w13_2 );
        delta4(i) = y14(i) * ( 1 - y14(i) ) * ( delta_out(i) * w14_2 );

        %Isejimo koeficientu atnaujinimas
        w11_2 = w11_2 + h * delta_out(i) * y11(i); 
        w12_2 = w12_2 + h * delta_out(i) * y12(i);
        w13_2 = w13_2 + h * delta_out(i) * y13(i);
        w14_2 = w14_2 + h * delta_out(i) * y14(i);

        b1_2 = b1_2 + h * delta_out(i);
    
        %Paslepto sluoksnio koeficientu atnaujinimas  
        w11 = w11 + h * delta1(i) * x(i);
        w21 = w21 + h * delta2(i) * x(i);
        w31 = w31 + h * delta3(i) * x(i);
        w41 = w41 + h * delta4(i) * x(i);
        w12 = w12 + h * delta1(i) * x(i);
        w22 = w22 + h * delta2(i) * x(i);
        w32 = w32 + h * delta3(i) * x(i);
        w42 = w42 + h * delta4(i) * x(i);
        
        b1(1) = b1(1) + h * delta1(i);
        b1(2) = b1(2) + h * delta2(i);
        b1(3) = b1(3) + h * delta3(i);
        b1(4) = b1(4) + h * delta4(i);
        
    end
    
for i = 1:1:length(x)
   
         v11(i) = x(i) * w11 + x(i) * w12 + b1(1);
         v12(i) = x(i) * w21 + x(i) * w22 + b1(2); 
         v13(i) = x(i) * w31 + x(i) * w32 + b1(3); 
         v14(i) = x(i) * w41 + x(i) * w42 + b1(4); 
    
    y11(i) = 1/(1+exp(-v11(i)));
    y12(i) = 1/(1+exp(-v12(i)));
    y13(i) = 1/(1+exp(-v13(i)));
    y14(i) = 1/(1+exp(-v14(i)));
    
    v21(i) = (y11(i) * w11_2) + (y12(i) * w12_2) + (y13(i) * w13_2) + (y14(i) * w14_2) + b1_2;
    
    y(i) = v21(i);
    e(i) = d(i) - y(i);
    ERROR = ERROR + abs(e(i));
end

end

plot(x, d, 'kp', x, y, 'rh');
grid on;












