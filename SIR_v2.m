clear all; close all;

%on recupere les donnees
[N,T,Z,F,Hfull,mX0,PX0,Qw,Rv,X] = simulationDonnees(0);

% nombres de particules
nbParticules = 100;

v = [0.0025 0.0025 0.0025 0.0025 0.0025 0.0025];
Qw = diag(v)

seuil = nbParticules/3;
% Initialisation des poids des particules
wPoids = ones(1, nbParticules) .* (1 / nbParticules);

% Generation des particules au premier instant

% R = chol(A) produces an upper triangular matrix R from the diagonal 
% and upper triangle of matrix A, satisfying the equation R'*R=A. 
% The chol function assumes that A is (complex Hermitian) symmetric. 
% If it is not, chol uses the (complex conjugate) transpose of 
% the upper triangle as the lower triangle. Matrix A must be positive definite.

% B = repmat(A,r1,...,rN) specifies a list of scalars, r1,..,rN, that describes 
% how copies of A are arranged in each dimension. 
% When A has N dimensions, the size of B is size(A).*[r1...rN]. 
% For example, repmat([1 2; 3 4],2,3) returns a 4-by-6 matrix.

Particules = chol(PX0)'*randn(6, nbParticules) + repmat(mX0, 1, nbParticules);
ParticulesNew = zeros(6,nbParticules);

hold on;

% Pour chaque instant de l'echantillonnage
for k = 2:N
    % Pour chaque particule
 % Mise a  jour des poids
    amers=1;
    if(isnan(Z(1,k)) && isnan(Z(3,k)))%on ne voit ni l'amer 1 ni la 2
        %si amers==0, pas de mise a jour des poids car pas de mesure
        amers=0;
    elseif(isnan(Z(1,k)))%on ne voit pas l'amer 1
        H=Hfull(3:4, :);
        R=Rv(3:4, 3:4);
        ZCurrent=Z(3:4,k);
    elseif(isnan(Z(3,k)))%on ne voit pas l'amer 2
        H=Hfull(1:2, :);
        R=Rv(1:2, 1:2);
        ZCurrent=Z(1:2,k);
    else %on voit les deux amers
        H=Hfull;
        R=Rv;
        ZCurrent=Z(:,k);
    end
    if(amers == 1)
        for i=1:nbParticules
            % Propagation de la particule xk-1 en xk
            ParticulesNew(:,i) = chol(Qw)'*randn(6, 1)+(F*Particules(:,i));
            % Calcul du poids si au moins un amer est percu
            wPoids(:,i) = wPoids(:,i) * (1/sqrt(det(2*pi*R))) * exp(-(1/2) * (ZCurrent-H*ParticulesNew(:,i))' * inv(R) * (ZCurrent-H*ParticulesNew(:,i)));   
        end
    else
        for i=1:nbParticules
            % Propagation de la particule xk-1 en xk
            ParticulesNew(:,i) = chol(Qw)'*randn(6, 1)+(F*Particules(:,i));
        end
    end
    
    %Normalisation des poids
    wPoids(:,:) = wPoids ./ sum(wPoids);
    
    %Calcul de la moyenne a posteriori
    esperance = zeros(6, 1);
    for i = 1 : nbParticules
        esperance = esperance + wPoids(i)*ParticulesNew(:,i);
    end
    
    %Calcul de la covariance a posteriori
    covariance = zeros(6,6);
    for i = 1 : nbParticules
        covariance = covariance + wPoids(i)*(ParticulesNew(:,i)-esperance)*(ParticulesNew(:,i)-esperance)';
    end
    
    %Reechantillonage pour eviter la degenerescence
    Neff = 1/sum(wPoids.^2);
    if Neff<seuil
        'Resample'
        [ParticulesNew,wPoids] = Resample(nbParticules, ParticulesNew, wPoids);
    end  
    
    %Affichage
    axis([-7,7,-7,7]);

    %Affichage des particules (points)
    particulesPlots = [];
    for j=1:nbParticules
        particulesPlots = [particulesPlots, plot(ParticulesNew(1,j),ParticulesNew(2,j),'.r'),plot(ParticulesNew(3,j),ParticulesNew(4,j),'.g'),plot(ParticulesNew(5,j),ParticulesNew(6,j),'.b')];
    end
    
    %Affichage des estimes des moyennes a posteriori (+ noirs) 
    P1 = plot(esperance(1),esperance(2),'+k');
    P2 = plot(esperance(3),esperance(4),'+k');
    P3 = plot(esperance(5),esperance(6),'+k');    

    %Affichage des des ellipses de confiances
    E1 = ellipse(esperance(1:2), covariance(1:2,1:2), 'r');
    E2 = ellipse(esperance(3:4), covariance(3:4,3:4), 'g');
    E3 = ellipse(esperance(5:6), covariance(5:6,5:6), 'b');
    
    %Affichage du vecteur d'etat reel (carres noirs)
    P1r = plot(X(1,k),X(2,k),'sk');%robot
    P2r = plot(X(3,k),X(4,k),'sk');%amer1
    P3r = plot(X(5,k),X(6,k),'sk');%amer2    
    
    pause(0.2);
    
    if(k~=N)
        delete([E1, P1, E2, P2, E3, P3,P1r,P2r,P3r]);
        delete(particulesPlots);
    end
    k
    Particules=ParticulesNew;
end