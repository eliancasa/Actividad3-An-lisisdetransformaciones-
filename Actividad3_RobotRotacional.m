%Limpieza de pantalla
clear all
close all
clc

%Declaración de variables simbólicas
syms th1(t) th2(t) th3(t) th4(t) th5(t) th6(t) t
syms d1 d4 d5 a2 a3

%Configuración del robot, 0 rotacional, 1 prismática
RP=[0 0 0 0 0 0];

%Vector de coordenadas articulares
Q= [th1 th2 th3 th4 th5 th6];
disp('Coordenadas generalizadas');
pretty(Q);

%Vector de velocidades generalizadas
Qp= diff(Q,t);
disp('Velocidades generalizadas');
pretty(Qp);

%Número de grados de libertad
GDL= size(RP,2);
GDL_str= num2str(GDL);

%% DEFINICIÓN DE ROTACIONES BÁSICAS

Rz = @(th)[cos(th) -sin(th) 0;
           sin(th)  cos(th) 0;
           0        0       1];

Ry = @(th)[cos(th)  0 sin(th);
           0        1 0;
          -sin(th)  0 cos(th)];

Rx = @(th)[1 0 0;
           0 cos(th) -sin(th);
           0 sin(th)  cos(th)];

%% Junta 1 
P(:,:,1)= [0;0;d1];
R(:,:,1)= Rz(th1);

%% Junta 2  
P(:,:,2)= [-a2;0;0];
R(:,:,2)= eye(3);

%% Junta 3  
P(:,:,3)= [0;a3;0];
R(:,:,3)= Ry(th3)*Rz(-pi/2);

%% Junta 4  
P(:,:,4)= [0;0;d4];
R(:,:,4)= Rx(th4);

%% Junta 5 
P(:,:,5)= [0;d5;0];
R(:,:,5)= Rx(th5);

%% Junta 6  
P(:,:,6)= [0;0;0];
R(:,:,6)= Rx(th6);

%Vector de ceros
Vector_Zeros= zeros(1,3);

%Inicialización
A(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
T(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
PO(:,:,GDL)= P(:,:,GDL); 
RO(:,:,GDL)= R(:,:,GDL); 
RO_inv(:,:,GDL)= R(:,:,GDL); 

for i = 1:GDL
    i_str= num2str(i);

    %Locales
    disp(strcat('Matriz de Transformación local A', i_str));
    A(:,:,i)=simplify([R(:,:,i) P(:,:,i); Vector_Zeros 1]);
    pretty(A(:,:,i));

    %Globales
    try
       T(:,:,i)= T(:,:,i-1)*A(:,:,i);
    catch
       T(:,:,i)= A(:,:,i);
    end

    disp(strcat('Matriz de Transformación global T', i_str));
    T(:,:,i)= simplify(T(:,:,i));
    pretty(T(:,:,i))

    RO(:,:,i)= T(1:3,1:3,i);
    RO_inv(:,:,i)= transpose(RO(:,:,i));
    PO(:,:,i)= T(1:3,4,i);
end

%% JACOBIANO LINEAL DIFERENCIAL
disp('Jacobiano lineal obtenido de forma diferencial');

jv_d = simplify(jacobian(PO(:,:,GDL), Q));

pretty(jv_d);

%% JACOBIANO ANALÍTICO

Jv_a = sym(zeros(3,GDL));
Jw_a = sym(zeros(3,GDL));

for k= 1:GDL
    
    if RP(k)==0
        
        if k==1
            z = sym([0;0;1]);
            p = PO(:,:,GDL);
        else
            z = RO(:,3,k-1);
            p = PO(:,:,GDL)-PO(:,:,k-1);
        end
        
        Jv_a(:,k)= simplify(cross(z,p));
        Jw_a(:,k)= simplify(z);
        
    elseif RP(k)==1
        
        if k==1
            Jv_a(:,k)= sym([0;0;1]);
        else
            Jv_a(:,k)= RO(:,3,k-1);
        end
        
        Jw_a(:,k)= sym([0;0;0]);
    end
end    

disp('Jacobiano lineal obtenido de forma analítica');
pretty(Jv_a);

disp('Jacobiano angular obtenido de forma analítica');
pretty(Jw_a);

%% VELOCIDADES

disp('Velocidad lineal obtenida mediante el Jacobiano lineal');
V=simplify(Jv_a*Qp');
pretty(V);

disp('Velocidad angular obtenida mediante el Jacobiano angular');
W=simplify(Jw_a*Qp');
pretty(W);