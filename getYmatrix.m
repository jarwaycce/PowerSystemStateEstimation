%---------------------------------------------------------------------------%
                  % 子程序 “getYmatrix.m”作用为计算节点导纳矩阵             
                  % 入口参数：节点参数矩阵bus,支路参数矩阵branch               
                  % 返回参数：简化节点导纳矩阵Ya(不考虑对地电容和变压器非标准变比，直接用支路电抗倒数）
                  %          完全节点导纳矩阵Yr(已考虑对地电容和变压器非标准变比）
                  %          参考节点编号nodeRe
%---------------------------------------------------------------------------%
function [Ya,Yr,nodeRe] = getYmatrix(bus,branch)

    nb=size(bus,1);                            % size(A,n) n=1返回行数，n=2返回列数
    nl=size(branch,1);
    Ya=zeros(nb,nb);                           % 对导纳矩阵赋初值0
    Yr=zeros(nb,nb);                           % 对导纳矩阵赋初值0
    for k=1:nl   
        I=branch(k,1);                         % 节点i
        J=branch(k,2);                         % 节点j
        Yt=1/(branch(k,3)+1j*branch(k,4));     % 支路参数 Z=R+jX                         
        Ym=1j*branch(k,5)/2;                   % 支路对地 jb/2
        K=branch(k,9);                         % i,j间变比K
        Yta=1/(1j*branch(k,4));                % 不考虑对地电容和变压器非标准变比，直接用支路电抗倒数
%% 计算Ya矩阵        
        Ya(I,I)=Ya(I,I)+Yta;                    
        Ya(J,J)=Ya(J,J)+Yta;
        Ya(I,J)=Ya(I,J)-Yta;
        Ya(J,I)=Ya(I,J);
%% 计算Yr矩阵
        if K~=0                                  % 变压器支路i:j=K:1
            Yr(I,I)=Yr(I,I)+Yt/K+(1-K)*Yt/K/K-Ym;% 考虑变压器非标准变比
            Yr(J,J)=Yr(J,J)+Yt/K+(K-1)*Yt/K-Ym;  % 变压器对地支路为负值
            Yr(I,J)=Yr(I,J)-Yt/K;                % 互导纳为负
        else                                     % 普通线路
            Yr(I,I)=Yr(I,I)+Yt+Ym;
            Yr(J,J)=Yr(J,J)+Yt+Ym;
            Yr(I,J)=Yr(I,J)-Yt;
        end
        Yr(J,I)=Yr(I,J);                         % 不考虑移相器支路，互导纳相同
    end
 %% 找到参考节点  
    nodeRe=0;                                 
    for t=1:nb
        Ys=(bus(t,5)+1j*bus(t,6))/100;     % 节点接地支路 Ys=Gs+jBs
        Yr(t,t)=Yr(t,t)+Ys;
        if bus(t,2)==3
            nodeRe=bus(t,1);
        end
    end
end
 