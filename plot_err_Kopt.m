pmin = 7; pmax = 12;
Nvec = 2.^(pmin:pmax)' + 1;
% data for 3 examples:
% chfield = 'd' with fac = 10,
% chfield = 'y' with fac = 4, chfield = 'y' with fac = 10
Kopt = [7,9,12,14,18,24;
        7,9,12,16,22,28;
        9,11,14,18,23,28];
Norm_ErrMax_olim2D = zeros(length(Nvec),3);
Norm_ERMS_olim2D = zeros(length(Nvec),3);
CPU_olim2D = zeros(length(Nvec),3);
mark = ['^','o','*'];

for j = 1 : 3   % index for 3 examples
    if j == 1
        chfield = 'd'; fac = 10;
    else
        if j == 2
            chfield = 'y'; fac = 4;
        else
            chfield = 'y'; fac = 10;
        end
    end
    for p = pmin : pmax
        i = p - pmin + 1;
        N = Nvec(i);
        fname = sprintf('measurements/olim2D_chfield_%c_fac_%d_N_%d_testK.txt',chfield,fac,N);
        data = load(fname);
        Norm_ErrMax_olim2D(i,j) = data(Kopt(j,i),2);
        Norm_ERMS_olim2D(i,j) = data(Kopt(j,i),3);
        CPU_olim2D(i,j) = data(Kopt(j,i),4);
    end
end 

%% plot CPU time vs N for all examples
figure(1); hold on; grid on;
for j = 1 : 3
    plot(Nvec,CPU_olim2D(:,j),'Linewidth',1,'Marker',mark(j));
end
hold off;
set(gca,'YScale','log','XScale','log');
xlabel('N'); ylabel('CPU time');
xticks(Nvec);
legend('chfield = ''d'', fac = 10','chfield = ''y'', fac = 4','chfield = ''y'', fac = 10','Location','Best');
% title('CPU time vs N');

%% plot Normalized ErrMax vs N for all examples
figure(2); hold on; grid on;
for j = 1 : 3
    plot(Nvec,Norm_ErrMax_olim2D(:,j),'Linewidth',1,'Marker',mark(j));
end
hold off;
set(gca,'YScale','log','XScale','log');
xlabel('N'); ylabel('Normalized ErrMax');
xticks(Nvec);
legend('chfield = ''d'', fac = 10','chfield = ''y'', fac = 4','chfield = ''y'', fac = 10','Location','Best');
% title('Normalized ErrMax vs N');

%% plot Normalized ERMS vs N for all examples
figure(3); hold on; grid on;
for j = 1 : 3
    plot(Nvec,Norm_ERMS_olim2D(:,j),'Linewidth',1,'Marker',mark(j));
end
hold off;
set(gca,'YScale','log','XScale','log');
xlabel('N'); ylabel('Normalized ERMS');
xticks(Nvec);
legend('chfield = ''d'', fac = 10','chfield = ''y'', fac = 4','chfield = ''y'', fac = 10','Location','Best');
% title('Normalized ERMS vs N');

%% plot CPU time vs Normalized ErrMax for all examples
figure(4); hold on; grid on;
for j = 1 : 3
    plot(Norm_ErrMax_olim2D(:,j),CPU_olim2D(:,j),'Linewidth',1,'Marker',mark(j));
end
hold off;
set(gca,'YScale','log','XScale','log');
xlabel('Normalized ErrMax'); ylabel('CPU time');
legend('chfield = ''d'', fac = 10','chfield = ''y'', fac = 4','chfield = ''y'', fac = 10','Location','Best');
% title('CPU time vs Normalized ErrMax');

%% plot CPU time vs Normalized ERMS for all examples
figure(5); hold on; grid on;
for j = 1 : 3
    plot(Norm_ERMS_olim2D(:,j),CPU_olim2D(:,j),'Linewidth',1,'Marker',mark(j));
end
hold off;
set(gca,'YScale','log','XScale','log');
xlabel('Normalized ERMS'); ylabel('CPU time');
legend('chfield = ''d'', fac = 10','chfield = ''y'', fac = 4','chfield = ''y'', fac = 10','Location','Best');
% title('CPU time vs Normalized ERMS');