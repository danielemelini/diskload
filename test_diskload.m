clear

% Load LNs

fid = fopen( 'REF_6371_loading_love_numbers_0_40000.txt', 'r' );
data = textscan( fid, '%d %f %f %f', 'headerlines', 1 );

h_love = data{2};
l_love = data{3};
k_love = data{4};

% Set some constants

alpha = 0.1;                           % Disk radius (degrees)
theta = linspace(0,alpha*5,100);       % Range of colatitudes wrt disk center
Tw    = 1;                             % Disk height (equivalent water height, m)
nmin  = 0;                             % Minimum degree
nmax  = [100:100:40000];               % Range of maximum degrees
imass = 0;                             % choose imass,0 or 1 (uncompensated/compensated load)
%
if imass==1
    fprintf('invoking a globally compensated load (icomp=1)\n')
else
    fprintf('invoking an uncompensated load (icomp=0)\n')
end

% Compute the disc response for the maximum value of nmax

[U,V,G]= diskload(alpha,imass,theta,Tw,nmin,nmax(end),h_love,k_love,l_love);

%% FIG 1

figure(1); clf;

plot(  theta./alpha, U, 'b', 'LineWidth', 1.5 );  hold on;
plot(  theta./alpha, V, 'r', 'LineWidth', 1.5 );
plot(  theta./alpha, G, 'g', 'LineWidth', 1.5 );

xlabel( '\theta/\alpha', 'FontSize', 16 );
ylabel( 'mm', 'FontSize', 16 );

xlim([0 5]);
ylim([-2.5 1]);

grid on;
legend( 'U', 'V', 'G', 'Location', 'best' );
tit1=['Disk radius \alpha = ',sprintf('%4.2f^\\circ',alpha)];
tit2=['  Load = ',sprintf('%4.2f',alpha),' m w.e.'];
title([tit1 tit2])

%% Examine sensitivity to nmax

% Compute the disc response vs nmax at theta = K * alpha
K=1.5;
[U,V,G]= diskload(alpha,imass,K*alpha,Tw,nmin,nmax,h_love,k_love,l_love);

%% FIG 2

figure(2); clf;

semilogx(  nmax, U, 'b', 'LineWidth', 1.5 );  hold on;
semilogx(  nmax, V, 'r', 'LineWidth', 1.5 );
semilogx(  nmax, G, 'g', 'LineWidth', 1.5 );

xlabel( 'n_{max}', 'FontSize', 16 );
ylabel( 'mm', 'FontSize', 16 );

xlim( [100 40000] );
ylim( [-1 .3] );     % good for K=1.5 or 1.25
%ylim( [-0.6 .2] );   % good for K = 1.75 or 2

nROT = 360 / alpha;
semilogx( [nROT nROT], ylim, 'k-.', 'LineWidth', 1.1 ); 
semilogx( 2*[nROT nROT], ylim, 'k-.', 'LineWidth', 1.1 ); 
semilogx(xlim,U(end)*[1 1],'b--', 'LineWidth', 0.5 );
semilogx(xlim,V(end)*[1 1],'r--', 'LineWidth', 0.5 );
semilogx(xlim,G(end)*[1 1],'g--', 'LineWidth', 0.5 );
hold off;


grid on;
legend( 'U', 'V', 'G', 'Location', 'best' );
Ht=title(['Truncation error: Loading response computed',...
       ' at \theta = ',sprintf('%4.2f',K),...
       ' \alpha as a function of n_{max} ']);
set(Ht,'FontSize',14)