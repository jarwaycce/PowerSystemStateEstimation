%---------------------------------------------------------------------------%
                  % 子程序 “gethmatrix1.m”作用为根据计算状态量计算h参数             
                  % 入口参数：节点参数矩阵bus，支路参数矩阵branch，
                  %          量测参数矩阵mdata，节点导纳矩阵Y
                  %          节点电压幅值ampV，节点电压相角angV 
                  % 返回参数：量测函数矢量h
%---------------------------------------------------------------------------%
function h = gethmatrix1(bus,branch,mdata,Y,angV,ampV)
    nbus=size(bus,1);                   
    nmdata=size(mdata,1);
    h = zeros(nmdata,1);
    for n=1:nmdata
        type=mdata(n,1);
        I=mdata(n,6);
        J=mdata(n,7);
        switch(type)
            case 0  % 节点电压
                h(n,1)=ampV(I,1);
            case 2  % 节点注入无功
                for s=1:nbus
                    h(n,1) = h(n,1)+ampV(I,1)*ampV(s,1)*(real(Y(I,s))*sin(angV(I,1)-angV(s,1))-imag(Y(I,s))*cos(angV(I,1)-angV(s,1)));
                end
            case 4  % 首端无功
                K=branch(mdata(n,8),9);
                gb=1/(branch(mdata(n,8),3)+1j*branch(mdata(n,8),4));
                yc=branch(mdata(n,8),5)/2;
                if K==0 % 线路支路
                    h(n,1) = -ampV(I,1)^2*(imag(gb)+yc)-ampV(I,1)*ampV(J,1)*(real(gb)*sin(angV(I,1)-angV(J,1))-imag(gb)*cos(angV(I,1)-angV(J,1)));
                else  % 变压器支路
                    h(n,1) = -1/(K^2)*ampV(I,1)^2*imag(gb)+1/K*ampV(I,1)*ampV(J,1)*imag(gb)*cos(angV(I,1)-angV(J,1));
                end
            case -4 % 末端无功
                K=branch(mdata(n,8),9);
                gb=1/(branch(mdata(n,8),3)+1j*branch(mdata(n,8),4));
                yc=branch(mdata(n,8),5)/2;
                if K==0 % 线路支路
                    h(n,1) = -ampV(J,1)^2*(imag(gb)+yc)+ampV(I,1)*ampV(J,1)*(real(gb)*sin(angV(I,1)-angV(J,1))+imag(gb)*cos(angV(I,1)-angV(J,1)));
                else  % 变压器支路
                    h(n,1) = -ampV(J,1)^2*imag(gb)+1/K*ampV(I,1)*ampV(J,1)*imag(gb)*cos(angV(I,1)-angV(J,1));
                end
            case 1 % 节点注入有功
                for s=1:nbus
                    h(n,1) = h(n,1)+ampV(I,1)*ampV(s,1)*(real(Y(I,s))*cos(angV(I,1)-angV(s,1))+imag(Y(I,s))*sin(angV(I,1)-angV(s,1)));
                end
            case 3 % 首端有功
                K=branch(mdata(n,8),9);
                gb=1/(branch(mdata(n,8),3)+1j*branch(mdata(n,8),4));
                if K==0 % 线路支路
                    h(n,1) = ampV(I,1)^2*real(gb)-ampV(I,1)*ampV(J,1)*(real(gb)*cos(angV(I,1)-angV(J,1))+imag(gb)*sin(angV(I,1)-angV(J,1)));
                else  % 变压器支路
                    h(n,1) = -1/K*ampV(I,1)*ampV(J,1)*imag(gb)*sin(angV(I,1)-angV(J,1));
                end
            case -3 %末端有功
                K=branch(mdata(n,8),9);
                gb=1/(branch(mdata(n,8),3)+1j*branch(mdata(n,8),4));
                if K==0 % 线路支路
                    h(n,1) = ampV(J,1)^2*real(gb)+ampV(I,1)*ampV(J,1)*(-real(gb)*cos(angV(I,1)-angV(J,1))+imag(gb)*sin(angV(I,1)-angV(J,1)));
                else  % 变压器支路
                    h(n,1) = 1/K*ampV(I,1)*ampV(J,1)*imag(gb)*sin(angV(I,1)-angV(J,1));
                end
        end
    end 
end
