chfield = 'y'; % 'd' or 'y'
fac = 4; % 4 or 10

pmin = 7; pmax = 12;
Nvec = 2.^(pmin:pmax)' + 1;
K = (1:40)';
Norm_ErrMax_olim2D = zeros(length(K),pmax-pmin+1);
Norm_ERMS_olim2D = zeros(length(K),pmax-pmin+1);
mark = ['^','o','*','s','d','v'];

for p = pmin : pmax
    i = p - pmin + 1;
    N = Nvec(i);
    fname = sprintf('measurements/olim2D_chfield_%c_fac_%d_N_%d_testK.txt',chfield,fac,N);
    data = load(fname);
    Norm_ErrMax_olim2D(:,i) = data(:,2);
    Norm_ERMS_olim2D(:,i) = data(:,3);
end

%% plot Normalized ErrMax vs K for all N
figure(1); hold on; grid on;
for p = pmin : pmax
    i = p - pmin + 1;
    plot(K,Norm_ErrMax_olim2D(:,i),'Linewidth',1,'Marker',mark(i));
end
hold off;
set(gca,'YScale','log');
xlabel('K'); ylabel('Normalized ErrMax');
legend('N = 129','N = 257','N = 513','N = 1025','N = 2049','N = 4097','Location','Best');
% ftitle = sprintf('case = ''%c'', fac = %d, Normalized ErrMax vs K',chfield,fac);
% title(ftitle);

%% plot Normalized ERMS vs K for all N
figure(2); hold on; grid on;
for p = pmin : pmax
    i = p - pmin + 1;
    plot(K,Norm_ERMS_olim2D(:,i),'Linewidth',1,'Marker',mark(i));
end
hold off;
set(gca,'YScale','log');
xlabel('K'); ylabel('Normalized ERMS');
legend('N = 129','N = 257','N = 513','N = 1025','N = 2049','N = 4097','Location','Best');
% ftitle = sprintf('case = ''%c'', fac = %d, Normalized ERMS vs K',chfield,fac);
% title(ftitle);