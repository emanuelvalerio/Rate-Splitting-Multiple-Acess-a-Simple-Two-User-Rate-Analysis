% This code was developed by EMANUEL VALERIO PEREIRA

% Given parameters
clc;clear all; clf
rand('state',0);
randn('state',0);
Nt = 2; % number of transmiter antennas (MISO system)
P = 1000;% Maximum power disponible
theta = linspace(0.000001,pi,50); % Adjust the number of points as needed
gamadb = linspace(-20, 0, 50); % Adjust the number of points as needed;
matrixT = zeros(length(gamadb),length(theta));
rhoMatrix = zeros(length(gamadb),length(theta));
idx = zeros(length(gamadb),length(theta));

ralationRateMatrix = zeros(length(gamadb),length(theta));
t0 = 0:0.1:1;
t = linspace(0,1,1);
bestRate = zeros(1,length(t));
tfoundMatrix = zeros(1,length(t0));
rateMatrix = zeros(1,length(t0));

for i = 1:length(gamadb)
    for j = 1:length(theta)
         gama = (10.^(gamadb(i)/20));
         [rho,h1,h2] = calculateRhoBasedOnGammaAndTheta(gama,theta(j)); % Calculate ρ
          rho = abs(rho);
          Gamma = (1/rho)*((1/norm(h2)^2)-(1/norm(h1)^2));     
          for m = 1:length(t0)
              [t_found,rate] = testePowerAllocated(Nt,P,h1,h2,rho,Gamma,gama,t0(m));
              tfoundMatrix(m) = t_found;
              rateMatrix(m) = rate;
          end
          t_found = max(tfoundMatrix);
          t_optimal = max(t_found,0);
          matrixT(i,j) = t_optimal;
          [powerAlocated , regime,~,idx,capacity] = newPowerAllocated(Nt,P,h1,h2,rho,Gamma,gama,t_optimal);
          rhoMatrix(i,j) = abs(rho);
          capacityMatrix(i,j)  = capacity;
          ind(i,j) = idx;

    end
end


[rateMulticast,rateOMA,rateNOMA,rateSDMA,rateRSMA] = rateSchemes(Nt,gamadb,theta,matrixT,P);

for i = 1:length(gamadb)
    for j = 1:length(theta)
        aux = 100*((rateRSMA(i,j)-max(rateSDMA(i,j),rateNOMA(i,j)))/max(rateSDMA(i,j),rateNOMA(i,j)));
        
           ralationRateMatrix(i,j) = max(0,aux);
       
    end
end

gamaPossibles = [-20,-10,-3,0];
for i = 1:length(gamaPossibles)
    for j = 1:length(theta)
         gama = (10.^(gamaPossibles(i)/20));
         [rho,h1,h2] = calculateRhoBasedOnGammaAndTheta(gama,theta(j)); % Calculate ρ
          rho = abs(rho);
          Gamma = (1/rho)*((1/norm(h2)^2)-(1/norm(h1)^2));     
          for m = 1:length(t0)
              [tFound,rate] = testePowerAllocated(Nt,P,h1,h2,rho,Gamma,gama,t0(m));
              tfoundMatrix(m) = tFound;
              rateMatrix(m) = rate;
          end
          tFound = max(tfoundMatrix);
          tOptimal = max(tFound,0);
          matrixT2(i,j) = tOptimal;
    end
end

[rateMulticast,rateOMA,rateNOMA,rateSDMA,rateRSMA] = rateSchemes(Nt,gamaPossibles,theta,matrixT2,P);
countOMA = 0;
for i = 1:length(gamaPossibles)
    for j = 1:length(theta)
        if(rateOMA(i,j)>rateNOMA(i,j)&&rateOMA(i,j)>rateMulticast&&rateOMA(i,j)>rateSDMA(i,j)&&rateOMA>rateRSMA(i,j))
            countOMA = countOMA+1;
        end
    end
end

% % figure 7
% OMAminus20 = 100*((sum(regimes(1,:) == "OMA"))/length(theta));
% NOMAminus20 = 100*((sum(regimes(1,:) == "NOMA"))/length(theta));
% RSMAminus20 = 100*((sum(regimes(1,:) == "RSMA"))/length(theta));
% SDMAminus20 = 100*((sum(regimes(1,:) == "SDMA"))/length(theta));
% multiCastminus20 = 100*((sum(regimes(1,:) == "Multicast"))/length(theta));
% 
% OMAminus10 = 100*((sum(regimes(2,:) == "OMA"))/length(theta));
% NOMAminus10 = 100*((sum(regimes(2,:) == "NOMA"))/length(theta));
% RSMAminus10 = 100*((sum(regimes(2,:) == "RSMA"))/length(theta));
% SDMAminus10 = 100*((sum(regimes(2,:) == "SDMA"))/length(theta));
% multiCastminus10 = 100*((sum(regimes(2,:) == "Multicast"))/length(theta));
% 
% OMAminus3 = 100*((sum(regimes(3,:) == "OMA"))/length(theta));
% NOMAminus3 = 100*((sum(regimes(3,:) == "NOMA"))/length(theta));
% RSMAminus3 = 100*((sum(regimes(3,:) == "RSMA"))/length(theta));
% SDMAminus3 = 100*((sum(regimes(3,:) == "SDMA"))/length(theta));
% multiCastminus3 = 100*((sum(regimes(3,:) == "Multicast"))/length(theta));
% 
% OMAgama0 = 100*((sum(regimes(4,:) == "OMA"))/length(theta));
% NOMAgama0 = 100*((sum(regimes(4,:) == "NOMA"))/length(theta));
% RSMAgama0 = 100*((sum(regimes(4,:) == "RSMA"))/length(theta));
% SDMAgama0 = 100*((sum(regimes(4,:) == "SDMA"))/length(theta));
% multiCastgama0 = 100*((sum(regimes(4,:) == "Multicast"))/length(theta));

figure(1)
num_levels = 30;
% Crie um gráfico de cores usando pcolor
h = pcolor(rhoMatrix(1,:), gamadb, matrixT);
colormap('parula'); % Escolha um mapa de cores, por exemplo, 'jet'
shading interp; % Interpola as cores para um visual mais suave

% Defina os limites dos eixos
xlim([min(rhoMatrix(1,:)), max(rhoMatrix(1,:))]);
ylim([min(gamadb), max(gamadb)]);

% Adicione uma barra de cores
ch = colorbar;

% Personalize o gráfico (opcional)
xlabel('\rho');
ylabel('Channel Strength disparity \gamma_{dB}[dB]');
title('Optimum t');
valores_x = [0, 0.2, 0.4, 0.6, 0.8, 1];
valores_y = [-20, -15, -10, -5, 0];
% Aplique os valores ao eixo x
xticks(valores_x);
yticks(valores_y);
set(ch, 'YTick', [0,0.2,0.4,0.6,0.8, 1]);
set(get(ch, 'ylabel'), 'string', 'Z');
% Adicionando contornos pretos
hold on; % Para manter o gráfico atual
contour(rhoMatrix(1,:),gamadb,matrixT,num_levels, 'k'); % 'k' especifica a cor preta para as linhas de contorno
hold off; % Liberar o gráfico

figure(2)

% Crie um gráfico de cores usando pcolor
h = pcolor(rhoMatrix(1,:), gamadb, ind);
colormap('parula'); % Escolha um mapa de cores, por exemplo, 'jet'
shading flat;

% Defina os limites dos eixos
xlim([min(rhoMatrix(1,:)), max(rhoMatrix(1,:))]);
ylim([min(gamadb), max(gamadb)]);

% Adicione uma barra de cores
ch = colorbar;

% Personalize o gráfico (opcional)
xlabel('\rho');
ylabel('Channel Strength disparity \gamma_{dB}[dB]');
title('Regions of operation for RSMA,SDMA,NOMA,OMA and Multicast.');
valores_x = [0, 0.2, 0.4, 0.6, 0.8, 1];
valores_y = [-20, -15, -10, -5, 0];
% Aplique os valores ao eixo x
xticks(valores_x);
yticks(valores_y);

figure(3)
num_levels = 10;
% Crie um gráfico de cores usando pcolor
h = pcolor(rhoMatrix(1,:), gamadb, ralationRateMatrix);
colormap('parula'); % Escolha um mapa de cores, por exemplo, 'jet'
shading interp; % Interpola as cores para um visual mais suave

% Defina os limites dos eixos
xlim([min(rhoMatrix(1,:)), max(rhoMatrix(1,:))]);
ylim([min(gamadb), max(gamadb)]);

% Adicione uma barra de cores
ch = colorbar;

% Personalize o gráfico (opcional)
xlabel('\rho');
ylabel('Channel Strength disparity \gamma_{dB}[dB]');
title('Fig 05 - Relative sum-rate gain [%] of RSMA over dynamic switching between SDMA and NOMA.');
valores_x = [0, 0.2, 0.4, 0.6, 0.8, 1];
valores_y = [-20, -15, -10, -5, 0];
% Aplique os valores ao eixo x
xticks(valores_x);
yticks(valores_y);
range = 0:1:(10*log10(P));
set(ch, 'YTick', range);
set(get(ch, 'ylabel'), 'string', 'Z');
% Adicionando contornos pretos
hold on; % Para manter o gráfico atual
contour(rhoMatrix(1,:),gamadb,ralationRateMatrix,num_levels, 'k'); % 'k' especifica a cor preta para as linhas de contorno
hold off; % Liberar o gráfico

% figure(4)
% percentage1 = [multiCastminus20,SDMAminus20,OMAminus20,NOMAminus20,RSMAminus20];
% percentage2 = [multiCastminus10,SDMAminus10,OMAminus10,NOMAminus10,RSMAminus10];
% percentage3 = [multiCastminus3,SDMAminus3,OMAminus3,NOMAminus3,RSMAminus3];
% percentage4 = [multiCastgama0,SDMAgama0,OMAgama0,NOMAgama0,RSMAgama0];
% bar([percentage1;percentage2;percentage3;percentage4], 'stacked');
% gamaValues = ["-20","-10","-3","0"];
% % Personalize o gráfico
% title('Gráfico de Barras Empilhadas de Porcentagens');
% xlabel('Porcentagem (%)');
% % yticks(1:length(categorias));
% % yticklabels(categorias);
% legend('Multicast', 'SDMA','OMA','NOMA','RSMA'); % Legenda para distinguir os conjuntos
% valores_x = [-20,-10,-3,0];
% % Aplique os valores ao eixo x
% xticks(valores_x);
% % % % Inverter a ordem das categorias para que a primeira categoria fique na parte superior
% %  set(gca, 'YDir', 'reverse');








