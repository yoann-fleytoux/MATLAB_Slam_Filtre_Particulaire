function [N,T,Z,F,Hfull,mX0,PX0,Qw,Rv,X] = simulationDonnees(plot_p);
%SIMULATIONDONNEES
%  Simulation d'une expérience : déplacement du robot, recueil des mesures, etc.
%  Syntaxe d'appel : [N,T,Z,F,Hfull,mX0,PX0,Qw,Rv,X] = simulationDonnees(plot_p);
%
%  Entrée :
%  . plot_p : si égal à 1, produit un affichage animé, sinon n'affiche rien
%
%  Sorties :
%  . N : nombre d'échantillons temporels
%  . T : vecteur des instants d'échantillonnage (1xN)
%  . Z : réalisation du processus aléatoire de mesure (4xN, éventuellement avec composantes de type NaN)
%  -> et comme ceci est de la simulation, sont également accessibles
%     . F, Hfull : matrices du modèle (respectivement (6x6), (4x6))
%     . mX0 : espérance du vecteur d'état à l'instant 0 (6x1)
%     . PX0 : covariance du vecteur d'état à l'instant 0 (6x6)
%     . Qw  : covariance du bruit de dynamique (supposé stationnaire) (6x6)
%     . Rv  : covariance du bruit de mesure (supposé stationnaire) (4x4)
%     . X   : réalisation du processus aléatoire d'état (6xN)
%
%  Z et X admettent AUTANT DE COLONNES QUE D'INSTANTS

N = 50;
deltaT = 1;
T = [0:N-1]*deltaT;
w=pi/4;

mX0 = [2;0;-1;4;4;-1];
PX0 = diag([.1 .2 .2 .2 .2 .2]);
Qw = diag([.05^2 .05^2 (1e-10)^2 (1e-10)^2 (1e-10)^2 (1e-10)^2]);
Rv = diag([.1^2 .1^2 .1^2 .1^2]);

F = blkdiag([cos(w*deltaT) -sin(w*deltaT) ; sin(w*deltaT) cos(w*deltaT)], eye(2), eye(2));

H1=[-1 0 1 0 0 0; 0 -1 0 1 0 0];
H2=[-1 0 0 0 1 0 ; 0 -1 0 0 0 1];
H=[H1;H2]; Hfull=H;
X = nan*ones(6,N); 
Z = nan*ones(4,N) ;
Source1Vue=zeros(1,N);
Source2Vue=zeros(1,N);

% Instant 0
W = chol(Qw)'*randn(6,N);
V = chol(Rv)'*randn(4,N);
X(:,1) = mX0+chol(PX0)'*randn(6,1);
% X = randn(n) returns an n-by-n matrix of normally distributed random numbers.

% R = chol(A) produces an upper triangular matrix R from the diagonal
% and upper triangle of matrix A, satisfying the equation R'*R=A. 
% The chol function assumes that A is (complex Hermitian) symmetric. 
% If it is not, chol uses the (complex conjugate) transpose of the upper triangle as the lower triangle. 
% Matrix A must be positive definite.

% Instants 1:(N-1);
for k=2:N                              
  X(:,k) = F*X(:,k-1) + W(:,k-1);
  
  if ((X(3,k)-X(1,k))>0)&&((X(4,k)-X(2,k))>0)
     Source1Vue(k)=1;
  end
  if ((X(5,k)-X(1,k))>0)&&((X(6,k)-X(2,k))>0)
     Source2Vue(k)=1;
  end
  
  if (Source1Vue(k)==1)&&(Source2Vue(k)==1)
    Z(:,k) = H*X(:,k) + V(:,k);
  elseif (Source1Vue(k)==1)&&(Source2Vue(k)==0)
    Z(:,k) = [H1*X(:,k) + V(1:2,k);nan(2,1)];
  elseif (Source1Vue(k)==0)&&(Source2Vue(k)==1)
    Z(:,k) = [nan(2,1);H2*X(:,k) + V(3:4,k)];
  else
    Z(:,k)=[nan(4,1)];
  end
  
  Z;
  X;
end;

if plot_p==1
    
scrsz=get(0,'ScreenSize');
currentfigure=figure('Position',[1 1 scrsz(3) scrsz(4)]);
nfig=1; 
set(nfig,'PaperPositionMode','auto'); 

for k=1:N
    h=subplot(1,2,2);
    set(h,'position',[0.53 0.05 .45 .92]);
    plot(X(1,k),X(2,k),'o','markeredgecolor','none','markerfacecolor','b');
    hold on
    grid on
    text(-4,-4,'VERITE TERRAIN (ESPACE D ETAT)')
    plot(X(1,max(k-15,1):k),X(2,max(k-15,1):k),'--b');
    visib_domain=[X(1,k) 5 5 X(1,k); X(2,k) X(2,k) 5 5];
    patch(visib_domain(1,:),visib_domain(2,:),'g','edgecolor','none','facealpha',.2);
    text(X(3,k)+.3,X(4,k),'Amer 1')
    text(X(5,k)+.3,X(6,k),'Amer 2')
    text(2,4.5,'CHAMP DE PERCEPTION','color','g')
    if Source1Vue(k)==1
        color='g';
        plot([X(3,k) X(3,k)],[X(2,k) X(4,k)],'--g','linewidth',2);
        text(X(3,k),X(2,k)-.2,sprintf('%s%s%s%s','(u^1-u^R)_','{',int2str(k),'}')); 
        plot([X(3,k) X(1,k)],[X(4,k) X(4,k)],'--g','linewidth',2);
        text(X(1,k)-.8,X(4,k),sprintf('%s%s%s%s','(v^1-v^R)_','{',int2str(k),'}')); 
    else
        color='r';
    end
    plot(X(3,k),X(4,k),'s','markeredgecolor','k','markerfacecolor',color);
    if Source2Vue(k)==1
        color='g';
        plot([X(5,k) X(5,k)],[X(2,k) X(6,k)],'--g','linewidth',2);
        text(X(5,k),X(2,k)-.2,sprintf('%s%s%s%s','(u^2-u^R)_','{',int2str(k),'}')); 
        plot([X(5,k) X(1,k)],[X(6,k) X(6,k)],'--g','linewidth',2);
        text(X(1,k)-.8,X(6,k),sprintf('%s%s%s%s','(v^2-v^R)_','{',int2str(k),'}')); 
    else
        color='r';
    end
    plot(X(5,k),X(6,k),'s','markeredgecolor','k','markerfacecolor',color);
    axis equal
    xlim([-5 5]);
    ylim([-5 5]);
    hold off
    
    g=subplot(1,2,1);
    set(g,'position',[0.03 0.05 .45 .92]);
    plot(Z(1,k),Z(2,k),'xk','markersize',7)
    hold on
    grid on
    text(5,5,'ESPACE DE MESURE')
    plot([0 0],[-.5 8],'k')
    plot([-.5 8],[0 0],'k')
    plot([0 Z(1,k)],[Z(2,k) Z(2,k)],'--k');
    plot([Z(1,k) Z(1,k)],[0 Z(2,k)],'--k');
    plot(Z(3,k),Z(4,k),'xk','markersize',7)
    plot([0 Z(3,k)],[Z(4,k) Z(4,k)],'--k');
    plot([Z(3,k) Z(3,k)],[0 Z(4,k)],'--k');
    if Source1Vue(k)==1
        text(Z(1,k),-.3,sprintf('%s%s%s%s','l^1_','{',int2str(k),'}'));
        text(-.5,Z(2,k),sprintf('%s%s%s%s','m^1_','{',int2str(k),'}')); 
    end
    if Source2Vue(k)==1
        text(Z(3,k),-.3,sprintf('%s%s%s%s','l^2_','{',int2str(k),'}'));
        text(-.5,Z(4,k),sprintf('%s%s%s%s','m^2_','{',int2str(k),'}')); 
    end
    axis equal
    xlim([-.5 8]);
    ylim([-.5 8]);
    hold off
    
    drawnow;
    pause(deltaT)
end

end
