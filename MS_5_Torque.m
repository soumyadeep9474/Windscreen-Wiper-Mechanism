clc; clear all;
L1=3.2;L2=20;L3=7.5;L0=5;L=20;
F_friction = 5;
t=0:0.05:7;
ang_speed =2;
theta = ang_speed*t;
A = [0;0];
D = L0*[-3;1];
F = D+[-5;0];
B = L1*[cos(theta);sin(theta)];
alpha = atan(-1/3)+pi;
theta2 = alpha-theta;
D_mag = 5*sqrt(10);
E = sqrt((D_mag*cos(alpha)-L1*cos(theta)).^2 + (D_mag*sin(alpha)- L1*sin(theta)).^2);
beta = acos((E.^2 + L2^2 - L3^2)./(2*E*L2));
gamma=-asin((D_mag*sin(alpha)-L1*sin(theta))./E)+pi;
omega=asin((L0-(L1*sin(theta)+L2*sin(gamma+beta)))/L3);
theta4 = pi-alpha+omega;
C = [(L1*cos(theta)+L2*cos(gamma+beta));(L1*sin(theta)+L2*sin(gamma+beta))];
delta = 0.874553;
E = D+[(L*cos(omega+delta));(L*sin(omega+delta))];
G = E+[-5;0];
E_up=D+[(L*cos(omega+delta));(L*sin(omega+delta))+10];
E_down=D+[(L*cos(omega+delta));(L*sin(omega+delta))-10];
E_x=E(1,:);
E_y=E(2,:);
E_vx = diff(E_x)./diff(t);
E_vy = diff(E_y)./diff(t);
E_v = sqrt(E_vx.^2 + E_vy.^2);
A_torque=(F_friction/ang_speed)*E_v;

Jacob= (L1/L3)*((L3*sin(theta2-theta4)+D_mag*sin(theta2)) ./ (L1*sin(theta2-theta4)+D_mag*sin(theta4)));
J_theta2 = (D_mag*L1*sin(theta4)/L3).*((D_mag*cos(theta2)+L3*cos(theta2-theta4)-L1)./((L1*sin(theta2-theta4)+D_mag*sin(theta4)).^2));
J_theta4 = (D_mag*L1*sin(theta2)/L3).*((D_mag*cos(theta4)+L1*cos(theta2-theta4)-L3)./((L1*sin(theta2-theta4)+D_mag*sin(theta4)).^2));
ang_acc_D=(ang_speed^2)*(J_theta2+Jacob.*J_theta4);

for i=1:length(t)
    ani = subplot(2,1,1);
    L1_bar=line([A(1) B(1,i)],[A(2) B(2,i)],'Color','green','linewidth',2);
    L2_bar=line([B(1,i) C(1,i)],[B(2,i) C(2,i)],'Color',[0.5 0 0.5 0.8],'linewidth',2);
    L3_bar=line([C(1,i) D(1)],[C(2,i) D(2)],'Color','c','linewidth',2);
    L_bar_fixed=line([E_down(1,i) E_up(1,i)],[E_down(2,i) E_up(2,i)],'Color',[0.8 1 0 0.25],'linewidth',3);
    L_bar=line([E_down(1,i) E_up(1,i)],[E_down(2,i) E_up(2,i)],'Color',[0.4940 0.1840 0.5560],'linewidth',2);
    L4_bar=line([D(1) E(1,i)],[D(2) E(2,i)],'Color','c','linewidth', 2);
    L5_bar=line([F(1) G(1,i)],[F(2) G(2,i)],'Color',[1 0 1],'linewidth',2);
    L6_bar=line([G(1,i) E(1,i)],[G(2,i) E(2,i)],'Color',[0.4940 0.1840 0.5560],'linewidth',2);
    A_circle = viscircles(A',0.15);
    B_circle = viscircles(B(:,i)',0.15);
    C_circle = viscircles(C(:,i)',0.15);
    D_circle = viscircles(D',0.15);
    E_circle = viscircles(E(:,i)',0.15);
    F_circle = viscircles(F',0.15);
    G_circle = viscircles(G(:,i)',0.15);
    axis(ani,'equal');
    set(gca,'XLim',[-35 5],'YLim',[-5 35]);
    str1 = 'E';
    str2 = ['Time elapsed' num2str(t(i)) 's'];
    time = text(-2,6,str2);
    pause(0.0000005);

    if(i < length(t))
        delete(A_circle);
        delete(B_circle);
        delete(C_circle);
        delete(D_circle);
        delete(E_circle);
        delete(F_circle);
        delete(G_circle);
        delete(L1_bar);
        delete(L2_bar);
        delete(L3_bar);
        delete(L4_bar);
        delete(L_bar);
        delete(L5_bar);
        delete(L6_bar);
        delete(time);
        vel = subplot(2,1,2);
        plot(vel,t(1:i),A_torque(1:i));
        set(vel,'XLim',[0 7],'YLim',[0 100]);
        xlabel(vel,'Time (s)');
        ylabel(vel,'Torque (N-m)');
        title(vel,'Input Torque at Crank A vs Time');
        grid on;

    end
end
