%Povilas_Kubilius_EKSfm-21
%Lab2
clc
clear
x = linspace(0.1,1,50);
d = ((1+ 0.6*sin(2/pi*x/0.7)) + 0.3*sin(2*pi*x))/2;
ERROR = 0;

%Atsitiktinis koeficientui generavimas

% Pirmasis sluoksnis

w1 = [randn(1) randn(1) randn(1) randn(1) randn(1)];
b1 = [randn(1) randn(1) randn(1) randn(1) randn(1)];

%Antrasis sluoksnis

w2 = [randn(1) randn(1) randn(1) randn(1) randn(1)];
b21= randn(1); 

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
    
    v11(i) = x(i) * w1(1) + b1(1);
    v12(i) = x(i) * w1(2) + b1(2); 
    v13(i) = x(i) * w1(3) + b1(3); 
    v14(i) = x(i) * w1(4) + b1(4); 
    v15(i) = x(i) * w1(5) + b1(5); 
    
    %Pirmojo sluoksnio isejimai  
    
    y11(i) = 1/(1+exp(-v11(i)));
    y12(i) = 1/(1+exp(-v12(i)));
    y13(i) = 1/(1+exp(-v13(i)));
    y14(i) = 1/(1+exp(-v14(i)));
    y15(i) = 1/(1+exp(-v15(i)));

    v21(i) = (y11(i) * w2(1))  +... 
    (y12(i) * w2(2)) + (y13(i) *... 
    w2(3)) + (y14(i) * w2(4))  + b21;

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

while ERROR > 0.05
    ERROR = 0;
    
    for i = 1:1:length(x)
        
        delta_out(i) = 1 * e(i);
        
        %Pasleptojo sluoksnio gradientas
        
        delta1(i) = y11(i) * ( 1 - y11(i) ) * ( delta_out(i) * w2(1) );
        delta2(i) = y12(i) * ( 1 - y12(i) ) * ( delta_out(i) * w2(2) );
        delta3(i) = y13(i) * ( 1 - y13(i) ) * ( delta_out(i) * w2(3) );
        delta4(i) = y14(i) * ( 1 - y14(i) ) * ( delta_out(i) * w2(4) );
        delta5(i) = y15(i) * ( 1 - y15(i) ) * ( delta_out(i) * w2(5) );

        %Isejimo koeficientu atnaujinimas
        
        w2(1) = w2(1) + h * delta_out(i) * y11(i); 
        w2(2) = w2(2) + h * delta_out(i) * y12(i);
        w2(3) = w2(3) + h * delta_out(i) * y13(i);
        w2(4) = w2(4) + h * delta_out(i) * y14(i);
        w2(5) = w2(5) + h * delta_out(i) * y15(i);
        b21 = b21 + h * delta_out(i);
    
        %Paslepto sluoksnio koeficientu atnaujinimas  
        
        w1(1) = w1(1) + h * delta1(i) * x(i);
        w1(2) = w1(2) + h * delta2(i) * x(i);
        w1(3) = w1(3) + h * delta3(i) * x(i);
        w1(4) = w1(4) + h * delta4(i) * x(i);
        w1(5) = w1(5) + h * delta5(i) * x(i);
        
        b1(1) = b1(1) + h * delta1(i);
        b1(2) = b1(2) + h * delta2(i);
        b1(3) = b1(3) + h * delta3(i);
        b1(5) = b1(5) + h * delta5(i);
        b1(4) = b1(4) + h * delta4(i);
        
    end
    
for i = 1:1:length(x)
   
    v11(i) = x(i) * w1(1) + b1(1);
    v12(i) = x(i) * w1(2) + b1(2); 
    v13(i) = x(i) * w1(3) + b1(3); 
    v14(i) = x(i) * w1(4) + b1(4); 
    v15(i) = x(i) * w1(5) + b1(5);
      
    y11(i) = 1/(1+exp(-v11(i)));
    y12(i) = 1/(1+exp(-v12(i)));
    y13(i) = 1/(1+exp(-v13(i)));
    y14(i) = 1/(1+exp(-v14(i)));
    y15(i) = 1/(1+exp(-v15(i)));
    
    v21(i) = (y11(i) * w2(1))  +... 
    (y12(i) * w2(2)) + (y13(i) *...
    w2(3)) + (y14(i) * w2(4)) + b21;
    
    y(i) = v21(i);
    e(i) = d(i) - y(i);
    ERROR = ERROR + abs(e(i));
    
end

end

plot(x, d, 'k-', x, y, 'r-');
grid on;












