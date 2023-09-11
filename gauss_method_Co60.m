function [b4, c4, b5, c5, b6, c6, plot_Co60] = gauss_method_Co60(file_isotope_Co60, file_fon)

f = fopen(file_fon,'r'); %окрываем файл с фоном
base = fread(f,inf,'single=>double'); %читаем в переменную массив значений
ambient_spectrum = base;
fclose(f); %закрываем файл

f = fopen(file_isotope_Co60,'r'); %окрываем файл со значениями изотопа 
isotope_spectrum1 = fread(f,inf, 'single=>double'); %читаем в переменную масив значений
fclose(f); %закрываем файл
len1 = min([size(ambient_spectrum, 1)/1024 size(isotope_spectrum1, 1)/1024]); % ищем минимальный по размеру массив из массивов изотопа и фона
ambient_spectrum = ambient_spectrum(1:len1*1024); %сокращем массив и приводим массивы к одной длине
isotope_spectrum1 = isotope_spectrum1(1:len1*1024); %сокращем массив и приводим массивы к одной длине
ambient_spectrum = sum (reshape(ambient_spectrum,1024,len1),2); %меняем размер массива фона
isotope_spectrum1 = sum (reshape(isotope_spectrum1,1024,len1),2); %меняем размер массива изотопа
isotope_spectrum1 = isotope_spectrum1-ambient_spectrum; %вычитаем фон из спектра изотопа 

plot_Co60 = isotope_spectrum1;
%figure;
%plot(isotope_spectrum1);
%grid on;

t = 1:1024;
t = t';

fitresult = createFit_Co60_1(t, isotope_spectrum1);
coeffs = coeffvalues(fitresult);
b4 = coeffs(2); c4 = coeffs(3);
b5 = coeffs(5); c5 = coeffs(6);

temp = b5+c5*2.335;
isotope_spectrum2 = isotope_spectrum1;
isotope_spectrum2(1:temp) = 0;

fitresult2 = createFit_Co60_2(t, isotope_spectrum2);
coeffs2 = coeffvalues(fitresult2);
b6 = coeffs2(2); c6 = coeffs2(3);
