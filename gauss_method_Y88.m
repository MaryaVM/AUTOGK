function [b7, c7, b8, c8, b9, c9, plot_Y88] = gauss_method_Y88(file_isotope_Y88, file_fon)

f = fopen(file_fon,'r'); %�������� ���� � �����
base = fread(f,inf,'single=>double'); %������ � ���������� ������ ��������
ambient_spectrum = base;
fclose(f); %��������� ����

f = fopen(file_isotope_Y88,'r'); %�������� ���� �� ���������� ������� 
isotope_spectrum1 = fread(f,inf, 'single=>double'); %������ � ���������� ����� ��������
fclose(f); %��������� ����
len1 = min([size(ambient_spectrum, 1)/1024 size(isotope_spectrum1, 1)/1024]); % ���� ����������� �� ������� ������ �� �������� ������� � ����
ambient_spectrum = ambient_spectrum(1:len1*1024); %�������� ������ � �������� ������� � ����� �����
isotope_spectrum1 = isotope_spectrum1(1:len1*1024); %�������� ������ � �������� ������� � ����� �����
ambient_spectrum = sum (reshape(ambient_spectrum,1024,len1),2); %������ ������ ������� ����
isotope_spectrum1 = sum (reshape(isotope_spectrum1,1024,len1),2); %������ ������ ������� �������
isotope_spectrum1 = isotope_spectrum1-ambient_spectrum; %�������� ��� �� ������� ������� 

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
