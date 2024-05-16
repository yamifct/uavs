% Para visualizar os gráficos, as tabelas têm de estar importadas no
% workspace. Para isso, basta arrastar os ficheiros cvs para o workspace a
% partir da diretoria onde estão inseridos e depois clicar em import em
% todos os ficheiros que aparecerem.
%% Analizar L2Data1
vars = who;
if ismember('L2Data1', vars)
    % Acelerômetro
    figure;
    subplot(3,1,1);
    plot(L2Data1.accx, 'r'); title('Acelerômetro - Eixo X');
    xlabel('Tempo'); ylabel('Acc X (m/s²)');

    subplot(3,1,2);
    plot(L2Data1.accy, 'g'); title('Acelerômetro - Eixo Y');
    xlabel('Tempo'); ylabel('Acc Y (m/s²)');

    subplot(3,1,3);
    plot(L2Data1.accz, 'b'); title('Acelerômetro - Eixo Z');
    xlabel('Tempo'); ylabel('Acc Z (m/s²)');

    % Plotar os dados do giroscópio
    figure;
    subplot(3,1,1);
    plot(L2Data1.gyrox, 'r'); title('Giroscópio - Eixo X');
    xlabel('Tempo'); ylabel('Gyro X (Graus/s)');

    subplot(3,1,2);
    plot(L2Data1.gyroy, 'g'); title('Giroscópio - Eixo Y');
    xlabel('Tempo'); ylabel('Gyro Y (Graus/s)');

    subplot(3,1,3);
    plot(L2Data1.gyroz, 'b'); title('Giroscópio - Eixo Z');
    xlabel('Tempo'); ylabel('Gyro Z (Graus/s)');

    % Calcular e exibir médias, desvios padrão e variância para o acelerômetro
    fprintf('Média do Acelerômetro - X: %f, Y: %f, Z: %f\n', mean(L2Data1.accx), mean(L2Data1.accy), mean(L2Data1.accz));
    fprintf('Desvio Padrão do Acelerômetro - X: %f, Y: %f, Z: %f\n', std(L2Data1.accx), std(L2Data1.accy), std(L2Data1.accz));
    fprintf('Covariância do Acelerômetro - X: %f, Y: %f, Z: %f\n', cov(L2Data1.accx), cov(L2Data1.accy), cov(L2Data1.accz));

    % Calcular e exibir médias , desvios padrão e variância para o giroscópio
    fprintf('Média do Giroscópio - X: %f, Y: %f, Z: %f\n', mean(L2Data1.gyrox), mean(L2Data1.gyroy), mean(L2Data1.gyroz));
    fprintf('Desvio Padrão do Giroscópio - X: %f, Y: %f, Z: %f\n', std(L2Data1.gyrox), std(L2Data1.gyroy), std(L2Data1.gyroz));
    fprintf('Covariância do Giroscópio - X: %f, Y: %f, Z: %f\n', cov(L2Data1.gyrox), cov(L2Data1.gyroy), cov(L2Data1.gyroz));

end

%% Analizar L2Data2

if exist('L2Data2', 'var')
    % Plotar os dados do acelerômetro
    figure;
    subplot(3,1,1);
    plot(L2Data2.accx, 'r'); title('Acelerômetro L2Data2 - Eixo X');
    xlabel('Tempo'); ylabel('Acc X (m/s²)');

    subplot(3,1,2);
    plot(L2Data2.accy, 'g'); title('Acelerômetro L2Data2 - Eixo Y');
    xlabel('Tempo'); ylabel('Acc Y (m/s²)');

    subplot(3,1,3);
    plot(L2Data2.accz, 'b'); title('Acelerômetro L2Data2 - Eixo Z');
    xlabel('Tempo'); ylabel('Acc Z (m/s²)');

    % Plotar os dados do giroscópio
    figure;
    subplot(3,1,1);
    plot(L2Data2.gyrox, 'r'); title('Giroscópio L2Data2 - Eixo X');
    xlabel('Tempo'); ylabel('Gyro X (Graus/s)');

    subplot(3,1,2);
    plot(L2Data2.gyroy, 'g'); title('Giroscópio L2Data2 - Eixo Y');
    xlabel('Tempo'); ylabel('Gyro Y (Graus/s)');

    subplot(3,1,3);
    plot(L2Data2.gyroz, 'b'); title('Giroscópio L2Data2 - Eixo Z');
    xlabel('Tempo'); ylabel('Gyro Z (Graus/s)');

    % Calcular e exibir médias, desvios padrão e covariância para o acelerômetro
    fprintf('Média do Acelerômetro L2Data2 - X: %f, Y: %f, Z: %f\n', mean(L2Data2.accx), mean(L2Data2.accy), mean(L2Data2.accz));
    fprintf('Desvio Padrão do Acelerômetro L2Data2 - X: %f, Y: %f, Z: %f\n', std(L2Data2.accx), std(L2Data2.accy), std(L2Data2.accz));
    fprintf('Covariância do Acelerômetro L2Data2 - X: %f, Y: %f, Z: %f\n', cov(L2Data2.accx), cov(L2Data2.accy), cov(L2Data2.accz));

    % Calcular e exibir médias, desvios padrão e covariância para o giroscópio
    fprintf('Média do Giroscópio L2Data2 - X: %f, Y: %f, Z: %f\n', mean(L2Data2.gyrox), mean(L2Data2.gyroy), mean(L2Data2.gyroz));
    fprintf('Desvio Padrão do Giroscópio L2Data2 - X: %f, Y: %f, Z: %f\n', std(L2Data2.gyrox), std(L2Data2.gyroy), std(L2Data2.gyroz));
    fprintf('Covariância do Giroscópio L2Data2 - X: %f, Y: %f, Z: %f\n', cov(L2Data2.gyrox), cov(L2Data2.gyroy), cov(L2Data2.gyroz));
    

pitch = asind(-L2Data2.accx/sqrt(L2Data2.accy.^2 + L2Data2.accz.^2+L2Data2.accx.^2));
roll = atan2(L2Data2.accy, L2Data2.accz);

figure;
subplot(2,1,1);
plot(pitch);
title('Pitch Angle');
ylabel('Graus');
xlabel('Tempo');

subplot(2,1,2);
plot(roll);
title('Roll Angle');
ylabel('Graus');
xlabel('Tempo');

end

%% Analizar L2Data3 e 4

if exist('L2Data3', 'var') && exist('L2Data4', 'var')
    % Cálculo dos ângulos de pitch e roll com base nos dados do acelerômetro para L2Data3
    pitch_L2Data3_calculated = asind(-L2Data3.accx/sqrt(L2Data3.accy.^2 + L2Data3.accz.^2+L2Data3.accx.^2));
    roll_L2Data3_calculated = atan2d(L2Data3.accy, L2Data3.accz);

    % Extração dos ângulos medidos diretamente pelo drone para L2Data3
    pitch_L2Data3_measured = L2Data3.stateEstimatepitch; 
    roll_L2Data3_measured = L2Data3.stateEstimateroll;

    % Cálculo dos ângulos de pitch e roll com base nos dados do acelerômetro para L2Data4
    pitch_L2Data4_calculated = asind(-L2Data4.accx/sqrt(L2Data4.accy.^2 + L2Data4.accz.^2+L2Data4.accx.^2));
    roll_L2Data4_calculated = atan2d(L2Data4.accy, L2Data4.accz);

    % Extração dos ângulos medidos diretamente pelo drone para L2Data4
    pitch_L2Data4_measured = L2Data4.stateEstimatepitch; 
    roll_L2Data4_measured = L2Data4.stateEstimateroll;

    % Plot comparativo para Pitch em L2Data3
    figure;
    subplot(2,1,1);
    plot(pitch_L2Data3_calculated, 'r'); hold on;
    plot(pitch_L2Data3_measured, 'b');
    title('Comparação de Pitch para L2Data3');
    legend('Calculado', 'Medido');
    ylabel('Pitch (graus)');
    xlabel('Tempo');

    subplot(2,1,2);
    plot(roll_L2Data3_calculated, 'r'); hold on;
    plot(roll_L2Data3_measured, 'b');
    title('Comparação de Roll para L2Data3');
    legend('Calculado', 'Medido');
    ylabel('Roll (graus)');
    xlabel('Tempo');

    % Plot comparativo para Pitch em L2Data4
    figure;
    subplot(2,1,1);
    plot(pitch_L2Data4_calculated, 'r'); hold on;
    plot(pitch_L2Data4_measured, 'b');
    title('Comparação de Pitch para L2Data4');
    legend('Calculado', 'Medido');
    ylabel('Pitch (graus)');
    xlabel('Tempo');

    subplot(2,1,2);
    plot(roll_L2Data4_calculated, 'r'); hold on;
    plot(roll_L2Data4_measured, 'b');
    title('Comparação de Roll para L2Data4');
    legend('Calculado', 'Medido');
    ylabel('Roll (graus)');
    xlabel('Tempo');

end

