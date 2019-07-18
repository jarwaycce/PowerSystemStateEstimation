%---------------------------------------------------------------------------%
                  % 子程序 “getBmatrix.m”作用为计算常数雅克比矩阵             
                  % 入口参数：简化节点导纳矩阵Ya，完全节点导纳矩阵Yr
                  %           支路参数矩阵branch，量测参数矩阵mdata，
                  %           参考节点nodeRe，参考节点电压初值ampV0
                  % 返回参数：P-Theta常数雅克比矩阵Ba
                  %          Q-V    常数雅克比矩阵Br  
%---------------------------------------------------------------------------%
function [Ba,Br] = getBmatrix(Ya,Yr,branch,mdata,nodeRe,ampV0)
    
    nbus=size(Ya,1);  
    nmdata=size(mdata,1);
    ImYa=imag(Ya);                           % 取Ya矩阵虚部
    ImYr=imag(Yr);                           % 取Yr矩阵虚部                    
    type=mdata(:,1);
    Pcount=length(find(type==1))+length(find(type==3))+length(find(type==-3));
    Qcount=nmdata-Pcount-1;                % 计算无功量测矢量个数（去除参考节点电压）
    Ba=zeros(Pcount,nbus-1);               % 有功矢量雅克比Ba矩阵初始化 ma×na阶（去除参考节点电压） 
    Br=zeros(Qcount,nbus-1);               % 无功矢量雅克比Br矩阵初始化 mr×nr阶（去除参考节点电压）
    Pcount=0;
    Qcount=0;
    for n=1:nmdata
        type=mdata(n,1);
        I=mdata(n,6);
        J=mdata(n,7);
        switch(type)
            case 0                     % 节点电压 Vi(i=j)
                if I~=nodeRe
                    Qcount=Qcount+1;
                    Br(Qcount,I)=-1/ampV0;  % 1=-v0*(-1/v0)
                end
            case 2                     % 节点注入无功 Qi(i=j)
                Qcount=Qcount+1;
                for t=1:nbus
                    if t==I
                        Br(Qcount,t)=ImYr(I,I)+sum(ImYr(I,:)); % (-Bii*v0^2+Qi)/v0=-v0*(Bii+sigma(Bij))
                    else
                        Br(Qcount,t)=ImYr(I,t);  % -Bij*v0=-v0*(Bij)
                    end
                end 
            case 4                     % 支路首端无功 Qij
                Qcount=Qcount+1;
                br=imag(1/(branch(mdata(n,8),3)+1j*branch(mdata(n,8),4)));
                y2c=branch(mdata(n,8),5);
                k=branch(mdata(n,8),9);
                if k~=0
                    Br(Qcount,I)= 2*br/k/k-br/k;
                    Br(Qcount,J)= -br/k;
                else
                    Br(Qcount,I)= br+y2c;  %-ImYr(I,J); % -v0*(b+2yc)   % b与B不同，B是互导纳的虚部，互导纳取得时候，Yij=-yij;
                    Br(Qcount,J)= -br; %ImYr(I,J);% v0*b=-v0*(-b)
                end
            case -4                    % 支路末端无功 Qji
                Qcount=Qcount+1;
                br=imag(1/(branch(mdata(n,8),3)+1j*branch(mdata(n,8),4)));
                y2c=branch(mdata(n,8),5);
                k=branch(mdata(n,8),9);
                if k~=0
                    Br(Qcount,I)= -br/k;
                    Br(Qcount,J)= 2*br-br/k;
                else
                    Br(Qcount,I)=-br;% ImYr(I,J); % v0*b=-v0*(-b)
                    Br(Qcount,J)=br+y2c; % -ImYr(I,J); % -v0*(b+2yc)
                end
            case 1                     % 节点注入有功Pi(i=j)
                Pcount=Pcount+1;
                for t=1:nbus
                    if t==I
                        Ba(Pcount,t)=ImYa(I,I)-sum(ImYa(I,:)); % -Bii*v0^2-Qi=-v0^2*(Bii-sigma(Bij))
                    else
                        Ba(Pcount,t)=ImYa(I,t); % -Bij*v0^2=-v0^2*(Bij)
                    end
                end        
            case 3                     % 支路首端有功Pij
                Pcount=Pcount+1;
                ba=-1/(branch(mdata(n,8),4));
                Ba(Pcount,I)= ba;% -ImYa(I,J); % -v0^2*b
                Ba(Pcount,J)=-ba;% ImYa(I,J);% v0^2*b
            case -3                   % 支路末端有功Pji
                Pcount=Pcount+1;
                ba=-1/(branch(mdata(n,8),4));
                Ba(Pcount,I)=-ba;% ImYa(I,J);% v0^2*b
                Ba(Pcount,J)=ba;% -ImYa(I,J);
        end
    end
    if nodeRe==1                         % 删去参考节点
        Br=Br(:,2:nbus);
        Ba=Ba(:,2:nbus);
    else
        Br=[Br(:,1:(nodeRe-1)) Br(:,((nodeRe+1):nbus))];
        Ba=[Ba(:,1:(nodeRe-1)) Ba(:,((nodeRe+1):nbus))];
    end
    
end
