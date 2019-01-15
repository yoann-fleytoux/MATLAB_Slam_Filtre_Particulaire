function [ ParticulesNew, WPoidsNew ] = Resample( nbParticules, Particules, WPoids )
    scp = cumsum(WPoids);
    ParticulesNew = zeros(size(Particules));
    WPoidsNew = zeros(size(WPoids));
    i=1;
    u1 = rand(1)*(1/nbParticules);
    u = zeros(nbParticules,1);
    for j=1:nbParticules
        u(j) = u1+(j-1)*(1/nbParticules);
        while u(j)>scp(i)
           i = i+1; 
        end
        ParticulesNew(:,j) = Particules(:,i);
        WPoidsNew(:,j) = 1/nbParticules;
    end
end

