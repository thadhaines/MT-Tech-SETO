% Generrate Figures from miniWECC test
% code functionalized

printFigs = 0;
% no Pss mac speed
compareMac_Spd( 'pstV2P3noPSSL.mat', 'pstV3p1noPSSL.mat', printFigs )
compareMac_Spd( 'pstV2P3noPSSL.mat', 'pstSETOnoPSSL.mat', printFigs )
compareMac_Spd( 'pstSETOnoPSSL.mat', 'pstV3p1noPSSL.mat', printFigs )

% pss mac speed
compareMac_Spd( 'pstV2P3PSSL.mat', 'pstV3p1PSSL.mat', printFigs )
compareMac_Spd( 'pstV2P3PSSL.mat', 'pstSETOPSSL.mat', printFigs )
compareMac_Spd( 'pstSETOPSSL.mat', 'pstV3p1PSSL.mat', printFigs )

% pss pss_out compare 
comparePSS_Out( 'pstV2P3PSSL.mat', 'pstV3p1PSSL.mat', printFigs )
comparePSS_Out( 'pstV2P3PSSL.mat', 'pstSETOPSSL.mat', printFigs )
comparePSS_Out( 'pstSETOPSSL.mat', 'pstV3p1PSSL.mat', printFigs )

%%
printFigs = 0;
compareBus_V( 'pstSETOnoPSSL.mat', 'pstV3p1noPSSL.mat', printFigs )
compareBus_Angle( 'pstSETOnoPSSL.mat', 'pstV3p1noPSSL.mat', printFigs )