function [graduirovka_graf, kalibrovka_graf, plot_fon, plot_Co57, plot_Co60, plot_Y88, plot_Cs137, plot_Eu152] = gauss_method_for_all_isotopes(file_fon, file_isotope_Co57, file_isotope_Co60, file_isotope_Y88, file_isotope_Cs137, file_isotope_Eu152)

%очищаем рабочее пространство
%clc;
%clear all;

%присваиваем переменным названия файлов для дальнейшего использования их в
%функциях
%file_fon = value_fon;
%file_isotope_Co57 = value_Co57;
%file_isotope_Co60 = value_Co60;
%file_isotope_Y88 = value_Y88;
%file_isotope_Cs137 = value_Cs137;
%file_isotope_Eu152 = value_Eu152;

file_fon_app = file_fon;
file_isotope_Co57_app = file_isotope_Co57;
file_isotope_Co60_app = file_isotope_Co60;
file_isotope_Y88_app = file_isotope_Y88;
file_isotope_Cs137_app = file_isotope_Cs137;
file_isotope_Eu152_app = file_isotope_Eu152;

%вызываем функции с аппроксимацией графиков в гауссианы и поиска
%коэффициентов с и b для каждого изотопа отдельно
[b1, c1, b2, c2, plot_fon] = gauss_method_K40_Tl208(file_fon_app);
[b3, c3, plot_Co57] = gauss_method_Co57(file_isotope_Co57_app, file_fon_app);
[b4, c4, b5, c5, b6, c6, plot_Co60] = gauss_method_Co60(file_isotope_Co60_app, file_fon_app);
[b7, c7, b8, c8, b9, c9, plot_Y88] = gauss_method_Y88(file_isotope_Y88_app, file_fon_app);
[b10, c10, plot_Cs137] = gauss_method_Cs137(file_isotope_Cs137_app, file_fon_app);
[b11, b12, c11, c12, plot_Eu152] = gauss_method_Eu152(file_fon_app, file_isotope_Eu152_app);

%создаем массивы энергий пиков, коэффициентов b и с соответственно
e = [1461, 2620, 122, 1173.2, 1332.5, 2505.7, 898, 1836, 2734, 661.6, 121.8, 344.3];
b_gauss = [b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12];
c_gauss = [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12];
c_gauss = c_gauss.*2.355; %умножаем на 2.355, чтобы получить 99% полуширины гауссианы, получится 1 сторона
c_gauss = c_gauss.*2; %умножаем на 2, чтобы получить ширину гауссианы полностью от левой до правой границы

%вызываем функцию для построения градуировочной шкалы
fitresult = graduirovka_gs(b_gauss, e);
coeffs = coeffvalues(fitresult); %достаем коэффициенты из fitresult
p1 = coeffs(1); %записываем коэффициенты в переменную 
p2 = coeffs(2); %записываем коэффициенты в переменную 

%записываем уравнение с 1024 каналами
s = 1024;
graduirovka_graf = zeros(s);
for x = 1:s
graduirovka_graf(x) = p1*x + p2;
end
graduirovka_graf = graduirovka_graf';
graduirovka_graf = sum(graduirovka_graf);

%строим график градуировочной шкалы
%figure();
%plot(graduirovka_graf);
%title('Градуировочная шкала');
%grid on;

fitresult2 = kalibrovka_gs(c_gauss, e);
coeffs2 = coeffvalues(fitresult2); %достаем коэффициенты из fitresult2
a = coeffs2(1); %записываем коэффициенты в переменную 
b = coeffs2(2); %записываем коэффициенты в переменную 

%записываем уравнение с 1024 каналами
kalibrovka_graf = zeros(s);
for x = 1:s
kalibrovka_graf(x) = a*x^b;
end
kalibrovka_graf = kalibrovka_graf';
kalibrovka_graf = sum(kalibrovka_graf);

%строим график калибровочной шкалы, аппрокимированной экспонентой 
%figure();
%plot(kalibrovka_graf);
%title('Калибровочная шкала');
%grid on;

%Ниже представлен алгоритм для построения калибровочной шкалы,
%аппроксимированной кривой мощности (power)

%вызываем функцию для построения калибровочной шкалы
%fitresult2 = kalibrovka(e, c_gauss);
%coeffs2 = coeffvalues(fitresult2); %достаем коэффициенты из fitresult2
%p1_2 = coeffs2(1); %записываем коэффициенты в переменную 
%p2_2 = coeffs2(2); %записываем коэффициенты в переменную 
%p3_2 = coeffs2(3); %записываем коэффициенты в переменную 

%записываем уравнение с 1024 каналами
%kalibrovka_graf = zeros(s);
%for x = 1:s
%kalibrovka_graf(x) = p1_2*x^2 + p2_2*x + p3_2;
%end
%kalibrovka_graf = kalibrovka_graf';
%kalibrovka_graf = sum(kalibrovka_graf);

%строим график калибровочной шкалы
%figure();
%plot(kalibrovka_graf);
%title('Калибровочная шкала');
%grid on;
