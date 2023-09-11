function [b1, c1, b2, c2, plot_fon] = gauss_method_K40_Tl208(file_fon)
% на вход (в круглых скобках) подается файл с записью фона
% на выход (к квадратных скобках) подается 4 коэффициента, по 2
% коэффициента (с и b) на каждый пик

f = fopen(file_fon,'r'); %окрываем файл с фоном
base = fread(f,inf,'single=>double'); %читаем в переменную массив значений
ambient_spectrum = base;
fclose(f); %закрываем файл

len1 = size(ambient_spectrum, 1)/1024; % ищем минимальный по размеру массив из массивов изотопа и фона
ambient_spectrum = ambient_spectrum(1:len1*1024); %сокращем массив и приводим массивы к одной длине
ambient_spectrum = sum (reshape(ambient_spectrum,1024,len1),2); %меняем размер массива фона

plot_fon = ambient_spectrum;
%строим график изотопа
%figure;
%plot(ambient_spectrum);
%grid on;

%создаем шкалу с 1024 каналами
t = 1:1024;
t = t'; %транспонируем массив, чтобы он был одного размера с массивом изотопа

fitresult = createFit_K40(t, ambient_spectrum);
coeffs = coeffvalues(fitresult);
b1 = coeffs(11);
c1 = coeffs(12);

temp = b1+c1*2.335; %вычисляем канал правой границы последнего найденного пика
ambient_spectrum2 = ambient_spectrum;
ambient_spectrum2(1:temp) = 0; %урезаем график от правой границы последнего найденного пика

fitresult2 = createFit_Tl208(t, ambient_spectrum2);
coeffs2 = coeffvalues(fitresult2);
b2 = coeffs2(2);
c2 = coeffs2(3);
