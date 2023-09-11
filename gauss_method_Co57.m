function [b3, c3, plot_Co57] = gauss_method_Co57(file_isotope_Co57, file_fon)
% �� ���� (� ������� �������) �������� ���� � ������� ���� � ���� � �������
% ������� Co57
% �� ����� (� ���������� �������) �������� 2 ������������(� � b)

f = fopen(file_fon,'r'); %�������� ���� � �����
base = fread(f,inf,'single=>double'); %������ � ���������� ������ ��������
ambient_spectrum = base;
fclose(f); %��������� ����

f = fopen(file_isotope_Co57,'r'); %�������� ���� �� ���������� ������� 
isotope_spectrum1 = fread(f,inf, 'single=>double'); %������ � ���������� ����� ��������
fclose(f); %��������� ����
len1 = min([size(ambient_spectrum, 1)/1024 size(isotope_spectrum1, 1)/1024]); % ���� ����������� �� ������� ������ �� �������� ������� � ����
ambient_spectrum = ambient_spectrum(1:len1*1024); %�������� ������ � �������� ������� � ����� �����
isotope_spectrum1 = isotope_spectrum1(1:len1*1024); %�������� ������ � �������� ������� � ����� �����
ambient_spectrum = sum (reshape(ambient_spectrum,1024,len1),2); %������ ������ ������� ����
isotope_spectrum1 = sum (reshape(isotope_spectrum1,1024,len1),2); %������ ������ ������� �������
isotope_spectrum1 = isotope_spectrum1-ambient_spectrum; %�������� ��� �� ������� ������� 

plot_Co57 = isotope_spectrum1;
%figure;
%plot(isotope_spectrum1);
%grid on;

t = 1:1024;
t = t';

fitresult3 = createFit_Co57(t, isotope_spectrum1);
coeffs = coeffvalues(fitresult3);
b3 = coeffs(2); c3 = coeffs(3);
