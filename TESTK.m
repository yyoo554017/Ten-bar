clear
clc
close

    r1=0.1;
    r2=0.05;
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
    
    T(1,1)=acosd(-1);
    T(2,1)=acosd(0);
    T(3,1)=acosd(1/sqrt(2));
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
     Kreduced=K(1:8,1:8);