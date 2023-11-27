clear;clc;
exact = @(x) sin(x);
exact_dx = @(x) cos(x);

f = @(x) sin(x);
g = sin(1);
h = -cos(0);

% number of element
N_el = [2 4 6 8 10 12 14 16];

% variable for storing error
e_L2 = zeros(length(N_el),1);
e_H1 = zeros(length(N_el),1);

for inel = 1 : length(N_el)
    %generate mesh
    n_el = N_el(inel);
    h_mesh = 1/n_el;
    x_coor = 0 : h_mesh/2 : 1;
    
    %IEN
    IEN = zeros(3,n_el);
    for e = 1 : n_el
        IEN(1,e) = 2*e-1;
        IEN(2,e) = 2*e;
        IEN(3,e) = 2*e+1;
    end
    
    %ID
    n_pt = 2*n_el+1;
    ID = zeros(1,n_pt);
    for i = 1 : n_pt
        ID(i) = i;
    end
    ID(n_pt) = 0;
    
    n_eq = n_pt - 1; % number of  equation
    
    % generate quadrature rule
    n_int = 6;
    [xi,weight] = Gauss(n_int,-1,1);
    K = zeros(n_eq,n_eq);
    F = zeros(n_eq,1);
    
    % Assembly of K and F
    for e = 1 : n_el
        
        k_e = zeros(3:3);
        f_e = zeros(3,1);
        
        x_ele = zeros(3,1);
        x_ele(1) = x_coor(IEN(1,e));
        x_ele(2) = x_coor(IEN(2,e));
        x_ele(3) = x_coor(IEN(3,e));
        
        for l = 1 : n_int
            dx_dxi = 0.0;
            x_l = 0.0;
            for a = 1 : 3
                dx_dxi = dx_dxi + x_ele(a) * QuadraticShape(a, xi(l),1);
                x_l = x_l +x_ele(a) * QuadraticShape(a,xi(l),0);
            end
            dxi_dx = 1/dx_dxi;
            
            for a = 1 : 3
                for b = 1 : 3
                    k_e(a,b) = k_e(a,b) + weight(l) * QuadraticShape(a,xi(l),1) * QuadraticShape(b,xi(l),1) * dxi_dx;
                end
            end
            
            for a = 1 : 3
                f_e(a) = f_e(a) + weight(l) * QuadraticShape(a,xi(l),0) * f(x_l) * dx_dxi;
            end 
        end
        
        % put element k and f into global K and F
        for a = 1 : 3
            AA = IEN(a,e);
            PP = ID(AA);
            if PP > 0
                F(PP) = F(PP)+f_e(a);
                for b = 1 : 3
                    BB = IEN(b,e);
                    QQ = ID(BB);
                    if QQ > 0
                        K(PP,QQ)=K(PP,QQ)+k_e(a,b);
                    else
                        F(PP) = F(PP) - k_e(a,b)*g;
                    end
                end
            end
        end
        if e == 1
            F(ID(IEN(1,e))) = F(ID(IEN(1,e)))+h;
        end
    end
    d = K \ F;
    disp = [d;g];
    
    %error
    e_l2 = 0.0;
    e_h1 = 0.0;
    for e = 1 : n_el
        for l = 1 : n_int
            uh = 0.0;
            uh_dx = 0.0;
            xl = 0.0;
            dx_dxi = 0.0;
            for a = 1 : 3
                dx_dxi = dx_dxi + x_coor(IEN(a,e)) * QuadraticShape(a, xi(l), 1);
            end
            
            for a = 1 : 3
                uh_dx = uh_dx + disp(IEN(a,e)) * QuadraticShape(a, xi(l), 1)/dx_dxi;
                uh = uh + disp(IEN(a,e)) * QuadraticShape(a, xi(l), 0);
                xl = xl + x_coor(IEN(a,e)) * QuadraticShape(a, xi(l), 0);
            end
            
            e_l2 = e_l2 + weight(l) * (uh - exact(xl))^2 * dx_dxi;
            e_h1 = e_h1 + weight(l) * (uh_dx - exact_dx(xl))^2 * dx_dxi;
        end
    end
    e_l2 = sqrt(e_l2) / sqrt(1/2-sin(2)/4);
    e_h1 = sqrt(e_h1) / sqrt(1/2+sin(2)/4);
    
    e_L2(inel) = e_l2;
    e_H1(inel) = e_h1;
end

%plot
hh = (1./N_el)';
plot(log(hh),log(abs(e_L2)));hold on;
plot(log(hh),log(abs(e_H1)));
xlabel('log(h)');
ylabel('log(error)');
title('Relative errors against the mesh size in log-log plot');
legend('eL2','eH1','Location','NorthWest');
%polyfit
t1 = polyfit(log(hh),log(e_L2),1);
t2 = polyfit(log(hh),log(e_H1),1);
fprintf('The scope of e_L2 is %.1f',t1(1));
fprintf('\n');
fprintf( 'The scope of e_H1 is %.1f',t2(1));
