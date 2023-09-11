function [b7, c7, b8, c8, b9, c9, plot_Y88] = gauss_method_Y88(file_isotope_Y88, file_fon)

f = fopen(file_fon,'r'); %окрываем файл с фоном
base = fread(f,inf,'single=>double'); %читаем в переменную массив значений
ambient_spectrum = base;
fclose(f); %закрываем файл

f = fopen(file_isotope_Y88,'r'); %окрываем файл со значениями изотопа 
isotope_spectrum1 = fread(f,inf, 'single=>double'); %читаем в переменную масив значений
fclose(f); %закрываем файл
len1 = min([size(ambient_spectrum, 1)/1024 size(isotope_spectrum1, 1)/1024]); % ищем минимальный по размеру массив из массивов изотопа и фона
ambient_spectrum = ambient_spectrum(1:len1*1024); %сокращем массив и приводим массивы к одной длине
isotope_spectrum1 = isotope_spectrum1(1:len1*1024); %сокращем массив и приводим массивы к одной длине
ambient_spectrum = sum (reshape(ambient_spectrum,1024,len1),2); %меняем размер массива фона
isotope_spectrum1 = sum (reshape(isotope_spectrum1,1024,len1),2); %меняем размер массива изотопа
isotope_spectrum1 = isotope_spectrum1-ambient_spectrum; %вычитаем фон из спектра изотопа 

plot_Y88 = isotope_spectrum1;
%figure;
%plot(isotope_spectrum1);
%grid on;

t = 1:1024;
t = t';

fitresult = createFit_Y88_1(t, isotope_spectrum1);
coeffs = coeffvalues(fitresult);
b7 = coeffs(2); c7 = coeffs(3);
b8 = coeffs(5); c8 = coeffs(6);

temp = b8+c8*2.335;
isotope_spectrum2 = isotope_spectrum1;
isotope_spectrum2(1:temp) = 0;

fitresult2 = createFit_Y88_2(t, isotope_spectrum2);
coeffs2 = coeffvalues(fitresult2);
b9 = coeffs2(2); c9 = coeffs2(3);
