% clear all
% close all
%  [x,y]=meshgrid(-100:1:-20, 3000:10:7000);
% z=exp(-7.4e-2 .*x+ 4.2e-3.* y - 12.6); 
% mesh(x,y,z)
%  xlabel('v_i'); ylabel('K_p'); zlabel('gKIR');
 
 clear all
close all
for k_p = [3000 6000 10000 13000 16000 20000];
 [v_i]=-110:0.1:-20;


 
  v_KIR_i = 4.5e-3 * k_p - 112;
  g_KIR_i = exp(-7.4e-2 .*v_i+ (4.2e-4).* k_p - 12.6);
  J_KIR_i =  7.5e2 * g_KIR_i / 1900 .* (v_i - v_KIR_i);
  plot(v_i,J_KIR_i)
 xlabel('v_i'); ylabel('J_KIR_i'); 
 legend('kp=3','kp= 6','kp= 10','kp= 13','kp= 16','kp= 20')
 grid on
 hold all
end