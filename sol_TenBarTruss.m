function [sigma, Q] = sol_TenBarTruss(r1, r2)
  
    % 定義各參數數值
    E=200*10^9;
    node=[18.28,9.14;18.28,0;9.14,9.14;9.14,0;0,9.14;0,0];
    A=zeros(10,1);
    A(1:6,1)=r1^2*pi();
    A(7:10,1)=r2^2*pi();
    L=zeros(10,1);
    L(1,1)=sqrt((node(5,1)-node(3,1))^2+(node(5,2)-node(3,2))^2);
    L(2,1)=sqrt((node(3,1)-node(1,1))^2+(node(3,2)-node(1,2))^2);
    L(3,1)=sqrt((node(6,1)-node(4,1))^2+(node(6,2)-node(4,2))^2); 
    L(4,1)=sqrt((node(4,1)-node(2,1))^2+(node(4,2)-node(2,2))^2);
    L(5,1)=sqrt((node(4,1)-node(3,1))^2+(node(4,2)-node(3,2))^2);
    L(6,1)=sqrt((node(2,1)-node(1,1))^2+(node(2,2)-node(1,2))^2);
    L(7,1)=sqrt((node(5,1)-node(4,1))^2+(node(5,2)-node(4,2))^2);
    L(8,1)=sqrt((node(6,1)-node(3,1))^2+(node(6,2)-node(3,2))^2);
    L(9,1)=sqrt((node(3,1)-node(2,1))^2+(node(3,2)-node(2,2))^2);
    L(10,1)=sqrt((node(4,1)-node(1,1))^2+(node(4,2)-node(1,2))^2);
    theta=zeros(10,1);
    theta(1,1)=acotd((node(5,1)-node(3,1))/(node(5,2)-node(3,2)));
    theta(2,1)=acotd((node(3,1)-node(1,1))/(node(3,2)-node(1,2)));
    theta(3,1)=acotd((node(6,1)-node(4,1))/(node(6,2)-node(4,2)));
    theta(4,1)=acotd((node(4,1)-node(2,1))/(node(4,2)-node(2,2)));
    theta(5,1)=acotd((node(4,1)-node(3,1))/(node(4,2)-node(3,2)));
    theta(6,1)=acotd((node(2,1)-node(1,1))/(node(2,2)-node(1,2)));
    theta(7,1)=acotd((node(5,1)-node(4,1))/(node(5,2)-node(4,2)));
    theta(8,1)=acotd((node(6,1)-node(3,1))/(node(6,2)-node(3,2)));
    theta(9,1)=acotd((node(3,1)-node(2,1))/(node(3,2)-node(2,2)));
    theta(10,1)=acotd((node(4,1)-node(1,1))/(node(4,2)-node(1,2)));
    % 開一個空白的剛性矩陣 (stiffness matrix)
    K=zeros(12);
     % 計算 stiffness matrix (可使用 add_element 函數)
    K = add_element(K,A(1,1),E,L(1,1),theta(1,1),3,5);
    K = add_element(K,A(2,1),E,L(2,1),theta(2,1),1,3);
    K = add_element(K,A(3,1),E,L(3,1),theta(3,1),4,6);
    K = add_element(K,A(4,1),E,L(4,1),theta(4,1),2,4);
    K = add_element(K,A(5,1),E,L(5,1),theta(5,1),3,4);
    K = add_element(K,A(6,1),E,L(6,1),theta(6,1),1,2);
    K = add_element(K,A(7,1),E,L(7,1),theta(7,1),4,5);
    K = add_element(K,A(8,1),E,L(8,1),theta(8,1),3,6);
    K = add_element(K,A(9,1),E,L(9,1),theta(9,1),2,3);
    K = add_element(K,A(10,1),E,L(10,1),theta(10,1),1,4);

    % 建立力矩陣
    F=[0 0 0 -1*10^7 0 0 0 -1*10^7 0 0 0 0]';
    Fr=F(1:8,1);
    % 建立空白位移矩陣
    Q =zeros(12,1);
  
    % 計算位移量 (F = KQ)
    Kr=K(1:8,1:8);
    %QF=inv(K)*F;
    Qr= inv(Kr)*Fr;
    Q(1:8,1)=Qr;
    % 建立空白應力矩陣
    sigma = zeros(10,1);
  
    % 計算應力 (stress) (可使用 compute_stress 函數)
    sigma(1,1)= compute_stress(Q,E,L(1,1),theta(1,1),3,5);
    sigma(2,1)= compute_stress(Q,E,L(2,1),theta(2,1),1,3);
    sigma(3,1)= compute_stress(Q,E,L(3,1),theta(3,1),4,6);
    sigma(4,1)= compute_stress(Q,E,L(4,1),theta(4,1),2,4);
    sigma(5,1)= compute_stress(Q,E,L(5,1),theta(5,1),3,4);
    sigma(6,1)= compute_stress(Q,E,L(6,1),theta(6,1),1,2);
    sigma(7,1)= compute_stress(Q,E,L(7,1),theta(7,1),4,5);
    sigma(8,1)= compute_stress(Q,E,L(8,1),theta(8,1),3,6);
    sigma(9,1)= compute_stress(Q,E,L(9,1),theta(9,1),2,3);
    sigma(10,1)= compute_stress(Q,E,L(10,1),theta(10,1),1,4);

    % (optional) compute reactions
    %KR=K(9:12,1:12);
    %R =KR*QF;
  
end