%---------------------------------------------------------------------------%
                  % 子程序 “ iteration.m”作用为迭代计算             
                  % 入口参数：节点参数矩阵bus,支路参数矩阵branch,
                  %          量测参数矩阵mdata,节点电压幅值、相角初值ampV0,angV0,
                  %          完全导纳矩阵Y,参考节点序号nodeRe
                  % 返回参数：节点电压幅值矢量ampV,节点电压相角矢量angV，
                  %          迭代次数iter，量测矢量z
                  %          量测函数矢量h
%---------------------------------------------------------------------------%
function [ampV,angV,iter,z,h] = iteration1(bus,branch,mdata,Y,ampV0,angV0,nodeRe)
    nbus=size(bus,1);
    nmdata=size(mdata,1);
    angV=zeros(nbus,1);                       
    ampV=zeros(nbus,1);       
    ampV(:,1)=ampV0;                          % 节点电压幅值赋初值v0
    angV(:,1)=angV0;                          % 节点电压相角赋初值0
    R=zeros(nmdata,nmdata);                   % 初始化权重矩阵
    z=zeros(nmdata,1);                        % 初始化量测量矩阵
    for i=1:nmdata                             
        R(i,i)=mdata(i,5);
        z(i,1)=mdata(i,2);
    end
    for n=1:100
        H = getJacmatrix1(branch,mdata,nodeRe,Y,ampV,angV);
        h = gethmatrix1(bus,branch,mdata,Y,angV,ampV);
        A = H'*R*H;
        B = H'*R*(z-h);
        del = A\B;
        if max(abs(del))<1e-5
            break;
        else
            ampV = ampV + del(1:nbus,1);
            temp = del(nbus+1:2*nbus-1,1);
            if nodeRe == 1
                angV(2:nbus,1) = angV(2:nbus,1) + temp;
            else
                angV(1:(nodeRe-1),1)    = angV(1:(nodeRe-1),1) + temp(1:(nodeRe-1),1);
                angV((nodeRe+1):nbus,1) = angV((nodeRe+1):nbus,1) + temp(nodeRe:nbus-1,1);
            end
        end
    end
    iter = n;
end