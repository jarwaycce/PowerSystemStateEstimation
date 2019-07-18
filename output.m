%---------------------------------------------------------------------------%
                  % 子程序 “ output.m”作用为输出界面显示             
                  % 入口参数：画图显示标志位draw_flag，状态量真值矩阵pfresult，
                  %          节点电压幅值矢量ampV,节点电压相角矢量angV，
                  %          迭代次数iter，迭代时间runtime，契合度ksi
                  % 返回参数：无
%---------------------------------------------------------------------------%
function []=output(alg_flag,draw_flag,pfresult,ampV,angV,iter,runtime,ksi)
    
    angV=180/pi*angV;
    time=datestr(now,'yyyy-mm-dd HH:MM');
    if alg_flag==1
        fprintf('                                       基于快速分解法的电力系统状态估计结果                                        \n');
        fprintf('                            Results of State Estimation Based on Fast Decoupled Method                           \n');
    else
        fprintf('                                   基于加权最小二乘法的电力系统状态估计结果                                        \n');
        fprintf('                       Results of State Estimation Based on Weight Least Squares Method                          \n');
    end
    fprintf('                                         %s    by J.W.Qi                                           \n',time);
    fprintf('                                             MATLAB Version：R2014a                                                  \n');
    fprintf('                                Intel(R) Core(TM) i7-6700HQ CPU @ 2.60GHz(8GB RAM)                                   \n');
    fprintf('==================================================================================================================\n');
    fprintf('节点号    节点电压幅值      节点电压幅值      幅值误差          节点电压相角       节点电压相角        相角误差   \n');
    fprintf('         (估计值/p.u.)     (真实值/p.u.)      (p.u.)         (估计值/degree)    (真实值/degree)      (degree)  \n');
    fprintf('==================================================================================================================\n');
    erampV=pfresult(:,2)-ampV(:,1);
    erangV=pfresult(:,3)-angV(:,1);
    for i=1:size(pfresult,1)
    fprintf('%6d %12.6f     %12.6f    %12.6f        %12.6f       %12.6f     %12.6f      \n',i,ampV(i,1),pfresult(i,2),erampV(i,1),angV(i,1),pfresult(i,3),erangV(i,1));
    end
    fprintf('==================================================================================================================\n');
    fprintf('估计效果评价指标：\n');
    fprintf('（1）算法的计算效率：【迭代次数】%-2d次                                 【运行时间】%-7.4f秒\n',iter,runtime);
    fprintf('（2）量测量估计精度：【契合程度】%4.2f%%\n',ksi);
    fprintf('（3）状态量估计精度: 【最大误差】电压幅值(p.u.)：%8.6f @ bus%4d    电压相角(degree)：%8.6f @ bus%4d\n',max(abs(erampV)),find(abs(erampV)==max(abs(erampV))),max(abs(erangV)),find(abs(erangV)==max(abs(erangV))));
    fprintf('                    【平均误差】电压幅值(p.u.)：%8.6f              电压相角(degree)：%8.6f\n',sum(abs(erampV))/size(pfresult,1),sum(abs(erangV))/size(pfresult,1));
    fprintf('------------------------------------------------------------------------------------------------------------------\n');
    if draw_flag==1   
        subplot(2,1,1);
        plot(1:size(ampV,1),ampV,'-r');
        grid on;      %添加网格
        hold on;
        plot(1:size(ampV,1),pfresult(:,2),':b');
        legend('估计值','真实值');
        axis([1,size(ampV,1),0.90,1.10]); 
    %     set(gca, 'XTick',[1:size(ampV,1)]);                  
        xlabel('节点序号');ylabel('节点电压幅值（p.u.）');
    %     title('基于快速分解法的电力系统状态估计结果（电压幅值）比较');

        subplot(2,1,2);
        plot(1:size(angV,1),angV,'-r');
        grid on;      %添加网格
        hold on;
        plot(1:size(angV,1),pfresult(:,3),':b');
        legend('估计值','真实值');
        axis([1 size(ampV,1) -90 90]);
    %     set(gca, 'XTick',[1:size(angV,1)]);                 
        xlabel('节点序号');ylabel('节点电压相角（degree）');
    %     title('基于快速分解法的电力系统状态估计结果（电压相角）比较');
    end
end