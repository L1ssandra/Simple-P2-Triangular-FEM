%hw5auto.m �Զ�����������

L = 4; % �ܹ��ļ������
Table = zeros(L,4); % ���

N0 = 5; % ��ʼ���ʷֶ���

%% ���
for m = 1:L
    
    N = N0*2^(m - 1); % ��������
    
    Nx = N; Ny = N;
    
    if N == 80
        draw_solution = 1;
    else
        draw_solution = 0;
    end
    
    tic
    main % ������
    toc
    
    Table(m,1) = L2; % ��L2�������
    Table(m,3) = H1; % ��H1�������
    
    if m > 1 % ȷ���о�����
        Table(m,2) = log(Table(m - 1,1)./Table(m,1))./log(2); % ����L2 Error�Ľ�
        Table(m,4) = log(Table(m - 1,3)./Table(m,3))./log(2); % ����H1 Error�Ľ�
    end
    
end