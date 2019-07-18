%---------------------------------------------------------------------------%
                  % 子程序 “gethmatrix.m”作用为计算量测函数h参数             
                  % 入口参数：节点参数矩阵bus，支路参数矩阵branch，
                  %           量测参数矩阵mdata，完全节点导纳矩阵Yr,
                  %           常数雅克比矩阵Ba,Br，节点电压相角矢量angV，
                  %           节点电压幅值矢量ampV，参考节点号nodeRe
                  %           有功无功量测函数计算标志位haflag
                  % 返回参数：有功量测函数矢量ha,无功量测函数矢量hr
%---------------------------------------------------------------------------%
function [ha,hr] = gethmatrix(bus,branch,mdata,Yr,Ba,Br,angV,ampV,nodeRe,haflag)
    nP=size(Ba,1);                    % 有功量测量个数
    nQ=size(Br,1);                    % 无功量测量个数（不含参考节点电压量）
    nmdata=size(mdata,1);
    nbus=size(bus,1);
    ha=zeros(nP,1);
    hr=zeros(nQ,1);
    Pcount=0;
    Qcount=0;
    if haflag==1                                      % 只计算有功部分
        for n=1:nmdata
            type=mdata(n,1);
            I=mdata(n,6);
            J=mdata(n,7);
            switch(type)
                case 1                                 % 节点注入有功
                    Pcount=Pcount+1;
                    for m=1:nbus
                        ha(Pcount,1) = ha(Pcount,1)+ampV(I,1)*ampV(m,1)*(real(Yr(I,m))*cos(angV(I,1)-angV(m,1))+imag(Yr(I,m))*sin(angV(I,1)-angV(m,1)));
                    end
                case 3                                 % 首端有功
                    Pcount=Pcount+1;
                    K=branch(mdata(n,8),9);
                    gb=1/(branch(mdata(n,8),3)+1j*branch(mdata(n,8),4));
                    if K==0                            % 线路支路
                        ha(Pcount,1) = ampV(I,1)^2*real(gb)-ampV(I,1)*ampV(J,1)*(real(gb)*cos(angV(I,1)-angV(J,1))+imag(gb)*sin(angV(I,1)-angV(J,1)));
                    else                               % 变压器支路
                        ha(Pcount,1) = -1/K*ampV(I,1)*ampV(J,1)*imag(gb)*sin(angV(I,1)-angV(J,1));
                    end
                case -3                                %末端有功
                    Pcount=Pcount+1;
                    K=branch(mdata(n,8),9);
                    gb=1/(branch(mdata(n,8),3)+1j*branch(mdata(n,8),4));
                    if K==0                            % 线路支路
                        ha(Pcount,1) = ampV(J,1)^2*real(gb)+ampV(I,1)*ampV(J,1)*(-real(gb)*cos(angV(I,1)-angV(J,1))+imag(gb)*sin(angV(I,1)-angV(J,1)));
                    else                               % 变压器支路
                        ha(Pcount,1) = 1/K*ampV(I,1)*ampV(J,1)*imag(gb)*sin(angV(I,1)-angV(J,1));
                    end
            end
        end
    else                                               % 只计算无功部分
        for n=1:nmdata
            type=mdata(n,1);
            I=mdata(n,6);
            J=mdata(n,7);
            switch(type)
                case 0                                 % 节点电压
                    if I~=nodeRe
                        Qcount=Qcount+1;
                        hr(Qcount,1)=ampV(I,1);
                    end
                case 2                                 % 节点注入无功
                    Qcount=Qcount+1;
                    for m=1:nbus
                        hr(Qcount,1) = hr(Qcount,1)+ampV(I,1)*ampV(m,1)*(real(Yr(I,m))*sin(angV(I,1)-angV(m,1))-imag(Yr(I,m))*cos(angV(I,1)-angV(m,1)));
                    end
                case 4                                 % 首端无功
                    Qcount=Qcount+1;
                    K=branch(mdata(n,8),9);
                    gb=1/(branch(mdata(n,8),3)+1j*branch(mdata(n,8),4));
                    yc=branch(mdata(n,8),5)/2;
                    if K==0                            % 线路支路
                        hr(Qcount,1) = -ampV(I,1)^2*(imag(gb)+yc)-ampV(I,1)*ampV(J,1)*(real(gb)*sin(angV(I,1)-angV(J,1))-imag(gb)*cos(angV(I,1)-angV(J,1)));
                    else                               % 变压器支路
                        hr(Qcount,1) = -1/(K^2)*ampV(I,1)^2*imag(gb)+1/K*ampV(I,1)*ampV(J,1)*imag(gb)*cos(angV(I,1)-angV(J,1));
                    end
                case -4                                % 末端无功
                    Qcount=Qcount+1;
                    K=branch(mdata(n,8),9);
                    gb=1/(branch(mdata(n,8),3)+1j*branch(mdata(n,8),4));
                    yc=branch(mdata(n,8),5)/2;
                    if K==0                            % 线路支路
                        hr(Qcount,1) = -ampV(J,1)^2*(imag(gb)+yc)+ampV(I,1)*ampV(J,1)*(real(gb)*sin(angV(I,1)-angV(J,1))+imag(gb)*cos(angV(I,1)-angV(J,1)));
                    else                               % 变压器支路
                        hr(Qcount,1) = -ampV(J,1)^2*imag(gb)+1/K*ampV(I,1)*ampV(J,1)*imag(gb)*cos(angV(I,1)-angV(J,1));
                    end
            end
        end
    end
end
