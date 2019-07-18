%---------------------------------------------------------------------------%
                  % 子程序 “getYmatrix1.m”作用为计算节点导纳矩阵             
                  % 入口参数：节点参数矩阵bus,支路参数矩阵branch               
                  % 返回参数：完全节点导纳矩阵Y(已考虑对地电容和变压器非标准变比）
                  %          参考节点编号nodeRe
%---------------------------------------------------------------------------%
function [Y,nodeRe] = getYmatrix1(bus,branch)

    nb=size(bus,1);                            % size(A,n) n=1返回行数，n=2返回列数
    nl=size(branch,1);
    Y=zeros(nb,nb);                           % 对导纳矩阵赋初值0
    for k=1:nl   
        I=branch(k,1);                         % 节点i
        J=branch(k,2);                         % 节点j
        Yt=1/(branch(k,3)+1j*branch(k,4));     % 支路参数 Z=R+jX                         
        Ym=1j*branch(k,5)/2;                   % 支路对地 jb/2
        K=branch(k,9);                         % i,j间变比K
%% 计算Yr矩阵
        if K~=0                                  % 变压器支路i:j=K:1
            Y(I,I)=Y(I,I)+Yt/K+(1-K)*Yt/K/K-Ym;% 考虑变压器非标准变比
            Y(J,J)=Y(J,J)+Yt/K+(K-1)*Yt/K-Ym;  % 变压器对地支路为负值
            Y(I,J)=Y(I,J)-Yt/K;                % 互导纳为负
        else                                     % 普通线路
            Y(I,I)=Y(I,I)+Yt+Ym;
            Y(J,J)=Y(J,J)+Yt+Ym;
            Y(I,J)=Y(I,J)-Yt;
        end
        Y(J,I)=Y(I,J);                         % 不考虑移相器支路，互导纳相同
    end
 %% 找到参考节点  
    nodeRe=0;                                 
    for t=1:nb
        Ys=(bus(t,5)+1j*bus(t,6))/100;     % 节点接地支路 Ys=Gs+jBs
        Y(t,t)=Y(t,t)+Ys;
        if bus(t,2)==3
            nodeRe=bus(t,1);
        end
    end
end
 