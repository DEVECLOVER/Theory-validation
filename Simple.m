clc
close all
clear all

SimpleFormular

function SimpleFormular()

    syms a1 a2 a3 a4 t1 t2 t3 f1 f2 c1 c2 u v d bt1 bt2 bt3
    
    base_form = (f2^2*cos(a2)^2*cos(a4)^2 + c2^2*cos(a4)^2*sin(a1)^2 + f1^2*cos(a1)^2*sin(a4)^2 + ...
        c1^2*sin(a1)^2*sin(a4)^2 + v^2*cos(a4)^2*sin(a1)^2 + u^2*sin(a1)^2*sin(a4)^2 - 2*c2*v*cos(a4)^2*sin(a1)^2 + ...
        c2^2*cos(a1)^2*cos(a4)^2*sin(a2)^2 - 2*c1*u*sin(a1)^2*sin(a4)^2 + c1^2*cos(a1)^2*sin(a2)^2*sin(a4)^2 + ...
        v^2*cos(a1)^2*cos(a4)^2*sin(a2)^2 + u^2*cos(a1)^2*sin(a2)^2*sin(a4)^2 + f1^2*sin(a1)^2*sin(a2)^2*sin(a4)^2 + ...
        2*c2*u*cos(a4)*sin(a1)^2*sin(a4) + 2*c1*v*cos(a4)*sin(a1)^2*sin(a4) + 2*f1*u*cos(a1)*sin(a1)*sin(a4)^2 - ...
        2*u*v*cos(a4)*sin(a1)^2*sin(a4) - 2*c2*v*cos(a1)^2*cos(a4)^2*sin(a2)^2 - ...
        2*c1*u*cos(a1)^2*sin(a2)^2*sin(a4)^2 - 2*c1*c2*cos(a4)*sin(a1)^2*sin(a4) - ...
        2*c1*f1*cos(a1)*sin(a1)*sin(a4)^2 - 2*c1*c2*cos(a1)^2*cos(a4)*sin(a2)^2*sin(a4) + ...
        2*c2*u*cos(a1)^2*cos(a4)*sin(a2)^2*sin(a4) + 2*c1*f1*cos(a1)*sin(a1)*sin(a2)^2*sin(a4)^2 + ...
        2*c1*v*cos(a1)^2*cos(a4)*sin(a2)^2*sin(a4) - 2*f1*u*cos(a1)*sin(a1)*sin(a2)^2*sin(a4)^2 - ...
        2*u*v*cos(a1)^2*cos(a4)*sin(a2)^2*sin(a4) + 2*c2*f1*cos(a1)*cos(a4)*sin(a1)*sin(a4) - ...
        2*f1*v*cos(a1)*cos(a4)*sin(a1)*sin(a4) + 2*c2*f2*cos(a1)*cos(a2)*cos(a4)^2*sin(a2) - ...
        2*f2*v*cos(a1)*cos(a2)*cos(a4)^2*sin(a2) - 2*c1*f2*cos(a1)*cos(a2)*cos(a4)*sin(a2)*sin(a4) + ...
        2*f2*u*cos(a1)*cos(a2)*cos(a4)*sin(a2)*sin(a4) - 2*f1*f2*cos(a2)*cos(a4)*sin(a1)*sin(a2)*sin(a4) - ...
        2*c2*f1*cos(a1)*cos(a4)*sin(a1)*sin(a2)^2*sin(a4) + 2*f1*v*cos(a1)*cos(a4)*sin(a1)*sin(a2)^2*sin(a4));

    midform = (c2-v)*cos(a4) - (c1-u)*sin(a4);
    F1 = f1*sin(a1)*sin(a2)*sin(a4);
    F2 = f2*cos(a2)*cos(a4);

    newform = (F2 - F1 + midform * sin(a2) * cos(a1))^2 + (f1*cos(a1)*sin(a4)+midform * sin(a1))^2;

    error_form = base_form - expand(newform)
    simplify(error_form)
end



function CalcNotConsistent()

%% 

    syms a1 a2 a3 t1 t2 t3 f1 f2 c1 c2 u v h w
    syms r11 r12 r13 r21 r22 r23 r31 r32 r33 t1 t2 t3 bt1 bt2 bt3
    Min = [f1 0 c1;0 f2 c2;0 0 1];  % camera internal params
    T_base = [r11 r12 r13 t1;r21 r22 r23 t2;r31 r32 r33 t3;0 0 0 1]; % external params
    
    Mout = [r11 r12 t1;r21 r22 t2;r31 r32 t3];

    Min*Mout


    R_x = [1 0 0 0;0 cos(a1) -sin(a1) 0;0 sin(a1) cos(a1) 0;0 0 0 1]; % X axis
    R_y = [cos(a2) 0 sin(a2) 0;0 1 0 0;-sin(a2) 0  cos(a2) 0;0 0 0 1]; % Y axis
    R_z = [cos(a3) -sin(a3) 0 0;sin(a3) cos(a3) 0 0;0 0 1 0;0 0 0 1]; % Z axis
    
    R_xy = (R_x * R_y) * R_z;  % eular angle

    R_xy(1:3,4) = [bt1 bt2 bt3];
%     R_xy = vpa(R_xy) 

    T_bias = T_base * R_xy;  
    
    R_base = T_base(1:3,[1,2,4]);  
    R_bias = T_bias(1:3,[1,2,4]);

    H_base = Min * R_base;   
    InvH_base = inv(H_base); 

    H_bias = Min * R_bias;  
    InvH_bias = inv(H_bias);

    
    Num = 2;
    baseworldpoint = [];
    biasworldpoint = [];
    U_lists = [0,h];
    V_lists = [0,w];
    for i = 1:2
        pp = [u + U_lists(i);v + V_lists(i);1];

        base_wp = InvH_base * pp;
        bias_wp = InvH_bias * pp;

        base_wp = base_wp / base_wp(3);
        bias_wp = bias_wp / bias_wp(3);
        
        base_wp = simplify(base_wp);
        bias_wp = simplify(bias_wp);

        baseworldpoint = [baseworldpoint,base_wp];
        biasworldpoint = [biasworldpoint,bias_wp];
    end
    basedis = baseworldpoint(:,end) - baseworldpoint(:,1);
    biasdis = biasworldpoint(:,end) - biasworldpoint(:,1);
    diffdis = biasdis - basedis  % difference error
    
%     partial_a3 = diff(bias(1),a3,1)
%     partial_bt3 = diff(bias(1),bt3,1)
    partial_bt1 = diff(diffdis(1),bt1,1);
    partial_bt2 = diff(diffdis(1),bt2,1);
  


%     R_x = R_base*R_x;
%     R_y = R_base*R_y;
%     R_z = R_base*R_z;
% 
%     R_x(:,3) = [t1;t2;t3];
%     R_y(:,3) = [t1;t2;t3];
%     R_z(:,3) = [t1;t2;t3];
% 
%     T_x = Min * R_x;
%     T_y = Min * R_y;
%     T_z = Min * R_z;
%     
%     T_xy = Min * R_x * R_y * R_z;
% 
%     Inv_Tx = inv(T_x); 
%     Inv_Ty = inv(T_y);
%     Inv_Tz = inv(T_z);
%     
%     Inv_Txy = inv(T_xy)
% 
%     Num = 6;
%     worldpoint = [];
%     for i = -Num:Num
%         pp = [u+i;v;1];
%         wp = Inv_Txy * pp;
%         wp = wp / wp(3);
%         wp = simplify(wp);
%         worldpoint = [worldpoint,wp];
%     end
%     
%     distances = worldpoint(:,2:end) - worldpoint(:,1:end-1);
%     
%     distances = simplify(distances);
%     
%     [N,D] = numden(distances);
% 
%     newdis = N ./ D;
%     diffdis = newdis - distances;
%     diffdis = simplify(diffdis);
% %     N = N ./ N(:,1);
%     simplify(N(1,1))
%     DD = simplify(D(1,1));
%     DD = expand(DD);
%     DD = simplify(DD)
%     D_diff = D(:,2:end) - D(:,1:end-1);
%     D_diff = simplify(D_diff);
%     
%     D_diff = simplify(D_diff - D_diff(:,1));
%     D_diff = simplify(D_diff)

    666
end

