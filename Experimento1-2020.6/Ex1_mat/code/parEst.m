function [sOut] = parEst(dW,dPrx,dPath)
  vtPrxmW = 10.^(dPrx/10);
  nSamples = length(dPrx);
  
  vtDesLargEst = [];
  vtDesPequeEst = [];
  %meia janela
  dMeiaJanela = round(dW/2);
  ij = 1;
  for ik = dMeiaJanela + 1: nSamples - dMeiaJanela
      vtDesLargEst(ij) = 10*log10(mean(vtPrxmW(ik-dMeiaJanela:ik+dMeiaJanela)));
      vtDesPequeEst(ij) = dPrx(ik) - vtDesLargEst(ij);
      ij = ij + 1;
  end
  
  %calculo da envoltória normalizada para posteriormente achar o fading
  in = dMeiaJanela +1 : nSamples - dMeiaJanela;
  vtPrxmWTemp = 10.^(dPrx(in)./10);
  desLargLin = 10.^(vtDesLargEst(1:length(in))./10);
  desLargLin = desLargLin';
  envNormal = sqrt(vtPrxmWTemp)./sqrt(desLargLin);
  

  %ajuste nos tamanhos devido a filtragem
  dPath = dPath(dMeiaJanela+1: nSamples-dMeiaJanela);
  vtPrx = dPrx(dMeiaJanela+1 : nSamples-dMeiaJanela);
  
  %calculo da reta de pathloss
  vtDistLog = log10(dPath);
  CoefReta = polyfit(vtDistLog, vtPrx, 1);
  NEst = -CoefReta(1)/10;
  %pathloss estimado
  vtPathLossEst = polyval(CoefReta, vtDistLog);
  %reshape para vetores linhas, já que foram criados automáticamente como vetores coluna no loop
  vtDesLargEst = vtDesLargEst';
  vtDesPequeEst = vtDesPequeEst';
  
  vtShadCorr = vtDesLargEst - vtPathLossEst;
  vtPathLossEst = - vtPathLossEst;
  %struct de retorno
  sOut.vtDesPequeEst = vtDesPequeEst;
  sOut.vtDesLargEst = vtDesLargEst;
  sOut.vtPathLossEst = vtPathLossEst;
  sOut.vtShadCorr = vtShadCorr;
  sOut.vtPrx = vtPrx;
  sOut.vtDistLog = vtDistLog;
  sOut.envNorm = envNormal;
  end
