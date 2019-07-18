%---------------------------------------------------------------------------%
                  % 子程序 “ iteration.m”作用为迭代计算             
                  % 入口参数：节点参数矩阵bus,支路参数矩阵branch,
                  %          量测参数矩阵mdata,节点电压幅值、相角初值ampV0,angV0,
                  %          完全导纳矩阵Yr,有功常数雅克比矩阵Ba,
                  %          无功常数雅克比矩阵Br，参考节点序号nodeRe
                  % 返回参数：节点电压幅值矢量ampV,节点电压相角矢量angV，
                  %          迭代次数iter，有功量测矢量Za,
                  %          无功量测矢量Zr,有功量测函数矢量ha,
                  %          无功量测函数矢量hr
%---------------------------------------------------------------------------%
function [ampV,angV,iter,Za,Zr,ha,hr] = iteration(bus,branch,mdata,ampV0,angV0,Yr,Ba,Br,nodeRe)
    mbus=size(bus,1);
    mP=size(Ba,1);
    mQ=size(Br,1);
    mr=size(mdata,1);
    
    angV=zeros(mbus,1);                    % 节点电压相角赋初值0
    ampV=zeros(mbus,1);       
    ampV(:,1)=ampV0;                       % 节点电压幅值赋初值v0
    angV(:,1)=angV0;  
    
    Ra=zeros(mP,mP);                       % 初始化权重矩阵
    Rr=zeros(mQ,mQ);
    Za=zeros(mP,1);                        % 初始化量测量矩阵
    Zr=zeros(mQ,1);
    Pcount=0;                              % 设置有功无功量计数标志
    Qcount=0;
    Pflag=0;                               % 置迭代标志
    Qflag=0;
    haflag=0;
%% 读入权重数据,量测数据  
    for i=1:mr                             
        type=mdata(i,1);
        if type==1||type==3||type==-3
            Pcount=Pcount+1;
            Ra(Pcount,Pcount)=mdata(i,5);
            Za(Pcount,1)=mdata(i,2);
        else
            if type==0 && mdata(i,6)==nodeRe && mdata(i,7)==nodeRe  % 参考节点不读
                continue;
            else
                Qcount=Qcount+1;
                Rr(Qcount,Qcount)=mdata(i,5);
                Zr(Qcount,1)=mdata(i,2);
            end
        end
    end
%% 计算常数信息矩阵
    A=ampV0^4*(-Ba)'*Ra*(-Ba);            
    B=ampV0^2*(-Br)'*Rr*(-Br);
%% 迭代过程
    for n=1:100
        if Pflag==0
%             ha=gethamatrix(bus,branch,mdata,Yr,Ba,angV,ampV);
            haflag=1;
            [ha,~] = gethmatrix(bus,branch,mdata,Yr,Ba,Br,angV,ampV,nodeRe,haflag);
            a=ampV0^2*(-Ba)'*Ra*(Za-ha);
            dAngV=A\a;            % 左除，A的逆乘以a
            if max(abs(dAngV))<1e-5
                Pflag=1;
            else                  % 修正相角
                if nodeRe==1
                    angV(2:mbus,1)         = angV(2:mbus,1)+dAngV(:,1);
                else
                    angV(1:(nodeRe-1),1)   = angV(1:(nodeRe-1),1)+dAngV(1:(nodeRe-1),1);
                    angV((nodeRe+1):mbus,1)= angV((nodeRe+1):mbus,1)+dAngV(nodeRe:(mbus-1),1);
                end
            end
        end
        if Qflag==0
%             hr= gethrmatrix(bus,branch,mdata,Yr,Br,angV,ampV,nodeRe);
            haflag=0;
            [~,hr] = gethmatrix(bus,branch,mdata,Yr,Ba,Br,angV,ampV,nodeRe,haflag);
            b=ampV0*(-Br)'*Rr*(Zr-hr);
            dAmpV=B\b;
            if max(abs(dAmpV))<(ampV0*1e-5)
                Qflag=1;
            else            % 修正幅值
                if nodeRe==1
                    ampV(2:mbus,1)          = ampV(2:mbus,1)+dAmpV(:,1);
                else 
                    ampV(1:(nodeRe-1),1)    = ampV(1:(nodeRe-1),1)+dAmpV(1:(nodeRe-1),1);
                    ampV((nodeRe+1):mbus,1) = ampV((nodeRe+1):mbus,1)+dAmpV(nodeRe:(mbus-1),1);
                end
            end
        end
        if Pflag==1 && Qflag==1
            break;
        end
    end
    iter=n;
end