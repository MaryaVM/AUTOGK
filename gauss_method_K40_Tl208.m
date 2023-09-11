function [b1, c1, b2, c2, plot_fon] = gauss_method_K40_Tl208(file_fon)
% �� ���� (� ������� �������) �������� ���� � ������� ����
% �� ����� (� ���������� �������) �������� 4 ������������, �� 2
% ������������ (� � b) �� ������ ���

f = fopen(file_fon,'r'); %�������� ���� � �����
base = fread(f,inf,'single=>double'); %������ � ���������� ������ ��������
ambient_spectrum = base;
fclose(f); %��������� ����

len1 = size(ambient_spectrum, 1)/1024; % ���� ����������� �� ������� ������ �� �������� ������� � ����
ambient_spectrum = ambient_spectrum(1:len1*1024); %�������� ������ � �������� ������� � ����� �����
ambient_spectrum = sum (reshape(ambient_spectrum,1024,len1),2); %������ ������ ������� ����

plot_fon = ambient_spectrum;
%������ ������ �������
%figure;
%plot(ambient_spectrum);
%grid on;

%������� ����� � 1024 ��������
t = 1:1024;
t = t'; %������������� ������, ����� �� ��� ������ ������� � �������� �������

fitresult = createFit_K40(t, ambient_spectrum);
coeffs = coeffvalues(fitresult);
b1 = coeffs(11);
c1 = coeffs(12);

temp = b1+c1*2.335; %��������� ����� ������ ������� ���������� ���������� ����
ambient_spectrum2 = ambient_spectrum;
ambient_spectrum2(1:temp) = 0; %������� ������ �� ������ ������� ���������� ���������� ����

fitresult2 = createFit_Tl208(t, ambient_spectrum2);
coeffs2 = coeffvalues(fitresult2);
b2 = coeffs2(2);
c2 = coeffs2(3);
