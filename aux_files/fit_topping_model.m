
a=2.71; % Lattice constant for Ru{0001}
n0 = (2/sqrt(3))*(1/a^2);
theta_exprimental = 0.023:0.001:0.031;

wf_data = f_wf_data; % Work Function data
indx = find(0.1 < wf_data(:,1) &  wf_data(:,1) < 0.5);

figure

x1=zeros(length(indx),2);
for i = 1:length(indx)

    subplot(3,ceil(length(indx)/3),i); hold on;
    plot(wf_data(1:indx(i),1),wf_data(1:indx(i),2),'ob');
    plot(wf_data(indx(i)+1:end,1),wf_data(indx(i)+1:end,2),'or');

    theta = wf_data(1:indx(i),1);
    wf_measured = wf_data(1:indx(i),2);

    wf_topping = @(x) dipole_dipole_repulsion.topping_work_function(x(1),n0,theta,x(2));

    x0 = [6 23];
    x1(i,:) = fminsearch(@(x)f(x,wf_topping,wf_measured),x0);
    plot(wf_data(1:indx(i),1),wf_topping(x1(i,:)))

    mu_eff = dipole_dipole_repulsion.effective_dipole_moment(x1(i,1),x1(i,2),n0,theta_exprimental);
    title(['\theta_{max}=' num2str(wf_data(indx(i)),2) ' , ' num2str(min(mu_eff),3) ' < \mu_{eff} <' num2str(max(mu_eff),3)])
end

function err = f(x,wf_topping,wf_measured)
if sum(x < 0) , err = 1e8; return; end

wf_calc = wf_topping(x);
err = sum((wf_calc-wf_measured).^2);
end

function wf_data = f_wf_data()
wf_data = [0.0057   -0.2422
    0.0057   -0.3674
    0.0158   -0.5010
    0.0158   -0.6430
    0.0302   -0.9269
    0.0374   -1.2192
    0.0503   -1.5449
    0.0776   -2.2046
    0.1020   -2.6305
    0.1681   -3.4154
    0.2112   -3.8163
    0.2615   -3.8163
    0.2859   -3.7745
    0.3204   -3.7077
    0.3376   -3.6994
    0.4052   -3.2985
    0.4856   -2.9729
    0.5129   -2.8727
    0.5546   -2.9144
    0.5690   -2.8977
    0.6293   -2.9979
    0.7414   -3.0313
    0.8032   -3.1148];
end