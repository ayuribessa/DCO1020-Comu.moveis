load('Prx_Real_2020_6.mat');

%%%%%%_______Início Parte 1 da Entrega_________%%%%%%

%distancia de meia janela é 5 para a primeira parte da entrega
dw = 5;
%assumindo Ptx = 1mw = 0dBm
txPower = 0;
[sOut] = parEst(dw,dPrx,dPath);
% plots
figure;
plot(sOut.vtDistLog, sOut.vtPrx); hold all;
plot(sOut.vtDistLog, txPower - sOut.vtPathLossEst,'linewidth',2)
plot(sOut.vtDistLog, txPower - sOut.vtPathLossEst + sOut.vtShadCorr,'LineWidth',2)
xlabel('log_{10}(d)');
ylabel('Potência [dBm]');
legend('Prx canal completo', 'Prx (somente perda de percurso)', 'Prx (perda de percurso + sombreamento)');
title('Prx original vs estimada');
% saveas(gcf,'results.png');

%%%%%%%_______Fim parte 1 da Entrega______%%%%%%


%%%%%%%_______Início Parte 2 da Entrega______%%%%%%

vtWs = [2 5 10];
for i=1:length(vtWs)
    sOut(i) = parEst(vtWs(i),dPrx,dPath);
    stdShad(i) = std(sOut(i).vtShadCorr);
    meanShad(i) = mean(sOut(i).vtShadCorr);
    coefReta = polyfit(sOut(i).vtDistLog, sOut(i).vtPrx, 1);
    nEst(i) = -coefReta(1)/10;
end
%Coloca em uma tabela para organização
varNames = {'W', 'stdShad', 'meanShad','nEst'};
T = table(vtWs',stdShad',meanShad',nEst','VariableNames',varNames);
disp('Média e desvio padrão do sombreamento para cada  valor de janela W:  ');
disp(T)

%%%%%%______Fim Entrega 02______%%%%%%


%%%%%%______Início Entrega 03______%%%%%%

sDistNames = [{'Weibull'};{'Rician'};{'Rayleigh'};{'Nakagami'}];
%Faz o ksteste armazena os valores de K e os parâmetros dos fits das
%dists
for i=1:length(vtWs)
%     disp(['    ------ Para janela W = ' num2str(vtWs(i)) ' ------']);
    for j = 1:length(sDistNames)
        data = [sOut(i).envNorm];
        pd(i,j)  = fitdist(data,sDistNames{j});
        x = linspace(min(data),max(data),length(data));
        tCDF = [x' cdf(pd(i,j),x)'];
        [h,p,k,c] = kstest(data,'CDF',tCDF);
        resulTest(j,i).h = h;
        resulTest(j,i).p = p;
        resulTest(j,i).k = k;
        resulTest(j,i).c = c;
        resulTest(j,i).name = sDistNames{j};
%         disp(['   Distribuição ' sDistNames{j} ': k = ' num2str(k) ', p-value = ' num2str(p)]);
        if (h==0)
            disp('Não rejeita a hipótese')
%         elseif(h==1)
%             disp('Rejeita a hipótese');
        end
    end
    disp(' ')
end

%Mostra a dist com menor K e seus parâmetros estimados pelo fitdist:
for i=1:length(vtWs)
    disp(['    ------------ Para W = ' num2str(vtWs(i)) ' ------------' ]);
%   [menorK,indice] = min([resulTest(:,i).k]);
    [out,indice] = sort([resulTest(:,i).k]);
    disp(['A melhor distribuição em função do menor K  pelo teste KS é: ' sDistNames{indice(1)} ', com k = ' num2str(out(1))]);
    disp(['E parâmetros: ' pd(i,indice(1)).ParameterNames ' iguais a: ']);
    disp(num2str(pd(i,indice(1)).ParameterValues));
    disp(' ');
    disp(['A segunda melhor distribuição em função do K é: ' sDistNames{indice(2)} ', com k = ' num2str(out(2))]);
    disp(['E parâmetros: ' pd(i,indice(2)).ParameterNames ' iguais a: ']);
    disp(num2str(pd(i,indice(2)).ParameterValues));
    disp(' ');
end
%%%%%%______ Fim Entrega 03______%%%%%%

%Achando parâmetros das dist. pelo pacote fitmethis   
% for i=1:length(vtWs)
%     disp(['--------Para W = ' num2str(vtWs(i)) '--------']);
%     fitwindow{i} = fitmethis([sOut(i).envNorm]','figure','off');
% end      

