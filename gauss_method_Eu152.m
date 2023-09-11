function [b11, b12, c11, c12, plot_Eu152] = gauss_method_Eu152(file_fon, file_isotope_Eu152)

f = fopen(file_fon,'r'); %окрываем файл с фоном
base = fread(f,inf,'single=>double'); %читаем в переменную массив значений
ambient_spectrum = base;
fclose(f); %закрываем файл

f = fopen(file_isotope_Eu152,'r'); %окрываем файл со значениями изотопа 
isotope_spectrum1 = fread(f,inf, 'single=>double'); %читаем в переменную масив значений
fclose(f); %закрываем файл
len1 = min([size(ambient_spectrum, 1)/1024 size(isotope_spectrum1, 1)/1024]); % ищем минимальный по размеру массив из массивов изотопа и фона
ambient_spectrum = ambient_spectrum(1:len1*1024); %сокращем массив и приводим массивы к одной длине
isotope_spectrum1 = isotope_spectrum1(1:len1*1024); %сокращем массив и приводим массивы к одной длине
ambient_spectrum = sum (reshape(ambient_spectrum,1024,len1),2); %меняем размер массива фона
isotope_spectrum1 = sum (reshape(isotope_spectrum1,1024,len1),2); %меняем размер массива изотопа
isotope_spectrum1 = isotope_spectrum1-ambient_spectrum; %вычитаем фон из спектра изотопа 

plot_Eu152 = isotope_spectrum1;
%figure;
%plot(isotope_spectrum1);
%grid on;

t = 1:1024;
t = t';

fitresult = createFit_Eu152_1_2(t, isotope_spectrum1);
coeffs1 = coeffvalues(fitresult);
b11 = coeffs1(5); c11 = coeffs1(6);
b12 = coeffs1(2); c12 = coeffs1(3);
