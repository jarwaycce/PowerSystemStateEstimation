%---------------------------------------------------------------------------%
                  % 子程序 “getJacmatrix1.m”作用为计算雅克比矩阵             
                  % 入口参数：支路参数矩阵branch，量测参数矩阵mdata
                  %           参考节点nodeRe，节点导纳矩阵Y，
                  %           节点电压幅值ampV，节点电压相角angV 
                  % 返回参数：雅克比矩阵H  
%---------------------------------------------------------------------------%
function H = getJacmatrix1(branch,mdata,nodeRe,Y,ampV,angV)
    
    nbus=size(Y,1);  
    nmdata=size(mdata,1);
    G=real(Y);
    B=imag(Y);                             % 取Y矩阵虚部
    H1=zeros(nmdata,nbus);                 % 关于状态量幅值ampV的偏导
    H2=zeros(nmdata,nbus);                 % 关于状态量相角angV的偏导
    for n=1:nmdata
        P=0;
        Q=0;
        type=mdata(n,1);
        I=mdata(n,6);
        J=mdata(n,7);
        switch(type)
            case 0                     % 节点电压 Vi(i=j)
                for m=1:nbus
                    if m == I
                        H1(n,m)=1;
                    else
                        H1(n,m)=0;
                    end
                   H2(n,m)=0;
                end
            case 1                     % 节点注入有功Pi(i=j)
                for m=1:nbus
                    P=P+ampV(I)*ampV(m)*(G(I,m)*cos(angV(I)-angV(m)) + B(I,m)*sin(angV(I)-angV(m)));
                    Q=Q+ampV(I)*ampV(m)*(G(I,m)*sin(angV(I)-angV(m)) - B(I,m)*cos(angV(I)-angV(m)));
                end
                for m=1:nbus
                    if m == I
                        H1(n,m)=(G(I,I)*ampV(I)^2+P)/ampV(I);
                        H2(n,m)=-B(I,I)*ampV(I)^2-Q;
                    else
                        H1(n,m)=ampV(I)*(G(I,m)*cos(angV(I)-angV(m)) + B(I,m)*sin(angV(I)-angV(m)));
                        H2(n,m)=ampV(I)*ampV(m)*(G(I,m)*sin(angV(I)-angV(m)) - B(I,m)*cos(angV(I)-angV(m)));
                    end
                end    
            case 2                     % 节点注入无功 Qi(i=j)
                for m=1:nbus
                    P=P+ampV(I)*ampV(m)*(G(I,m)*cos(angV(I)-angV(m)) + B(I,m)*sin(angV(I)-angV(m)));
                    Q=Q+ampV(I)*ampV(m)*(G(I,m)*sin(angV(I)-angV(m)) - B(I,m)*cos(angV(I)-angV(m)));
                end
                for m=1:nbus
                    if m == I
                        H1(n,m)=(-B(I,I)*ampV(I)^2+Q)/ampV(I);
                        H2(n,m)=-G(I,I)*ampV(I)^2+P;
                    else
                        H1(n,m)=ampV(I)*(G(I,m)*sin(angV(I)-angV(m)) - B(I,m)*cos(angV(I)-angV(m)));
                        H2(n,m)=-ampV(I)*ampV(m)*(G(I,m)*cos(angV(I)-angV(m)) + B(I,m)*sin(angV(I)-angV(m)));
                    end
                end 
            case 3                     % 支路首端有功Pij
                y = 1/(branch(mdata(n,8),3)+1j*(branch(mdata(n,8),4)));
                g = real(y);
                b = imag(y);
                k = branch(mdata(n,8),9);
                if k~=0
                    H1(n,I)=-(ampV(J)*b*sin(angV(I)-angV(J)))/k;
                    H1(n,J)=-(ampV(I)*b*sin(angV(I)-angV(J)))/k;
                    H2(n,I)=-(ampV(I)*ampV(J)*b*cos(angV(I)-angV(J)))/k;
                    H2(n,J)=(ampV(I)*ampV(J)*b*cos(angV(I)-angV(J)))/k;
                else
                    H1(n,I)=2*ampV(I)*g-ampV(J)*g*cos(angV(I)-angV(J))-ampV(J)*b*sin(angV(I)-angV(J));
                    H1(n,J)=-ampV(I)*(g*cos(angV(I)-angV(J))+b*sin(angV(I)-angV(J)));
                    H2(n,I)=ampV(I)*ampV(J)*(g*sin(angV(I)-angV(J))-b*cos(angV(I)-angV(J)));
                    H2(n,J)=-ampV(I)*ampV(J)*(g*sin(angV(I)-angV(J))-b*cos(angV(I)-angV(J)));
                end
            case 4                     % 支路首端无功 Qij
                y = 1/(branch(mdata(n,8),3)+1j*(branch(mdata(n,8),4)));
                yc = branch(mdata(n,8),5)/2;
                g = real(y);
                b = imag(y);
                k = branch(mdata(n,8),9);
                if k~=0
                    H1(n,I)=-(2*ampV(I)*b)/k/k + (ampV(J)*b*cos(angV(I)-angV(J)))/k;
                    H1(n,J)=(ampV(I)*b*cos(angV(I)-angV(J)))/k;
                    H2(n,I)=-(ampV(I)*ampV(J)*b*sin(angV(I)-angV(J)))/k;
                    H2(n,J)=(ampV(I)*ampV(J)*b*sin(angV(I)-angV(J)))/k;
                else
                    H1(n,I)=-2*ampV(I)*(b+yc)-ampV(J)*(g*sin(angV(I)-angV(J))-b*cos(angV(I)-angV(J)));
                    H1(n,J)=-ampV(I)*(g*sin(angV(I)-angV(J))-b*cos(angV(I)-angV(J)));
                    H2(n,I)=-ampV(I)*ampV(J)*(g*cos(angV(I)-angV(J))+b*sin(angV(I)-angV(J)));
                    H2(n,J)=ampV(I)*ampV(J)*(g*cos(angV(I)-angV(J))+b*sin(angV(I)-angV(J)));
                end
            case -3                    % 支路末端有功Pji
                y = 1/(branch(mdata(n,8),3)+1j*(branch(mdata(n,8),4)));
                g = real(y);
                b = imag(y);
                k = branch(mdata(n,8),9);
                if k~=0
                    H1(n,I)=(ampV(J)*b*sin(angV(I)-angV(J)))/k;
                    H1(n,J)=(ampV(I)*b*sin(angV(I)-angV(J)))/k;
                    H2(n,I)=(ampV(I)*ampV(J)*b*cos(angV(I)-angV(J)))/k;
                    H2(n,J)=-(ampV(I)*ampV(J)*b*cos(angV(I)-angV(J)))/k;
                else
                    H1(n,I)=ampV(J)*(-g*cos(angV(I)-angV(J))+b*sin(angV(I)-angV(J)));
                    H1(n,J)=2*ampV(J)*g+ampV(I)*(-g*cos(angV(I)-angV(J))+b*sin(angV(I)-angV(J)));
                    H2(n,I)=ampV(I)*ampV(J)*(g*sin(angV(I)-angV(J))+b*cos(angV(I)-angV(J)));
                    H2(n,J)=-ampV(I)*ampV(J)*(g*sin(angV(I)-angV(J))+b*cos(angV(I)-angV(J)));
                end
            case -4                    % 支路末端无功 Qji
                y = 1/(branch(mdata(n,8),3)+1j*(branch(mdata(n,8),4)));
                yc = branch(mdata(n,8),5)/2;
                g = real(y);
                b = imag(y);
                k = branch(mdata(n,8),9);
                if k~=0
                    H1(n,I)=(ampV(J)*b*cos(angV(I)-angV(J)))/k;
                    H1(n,J)=-2*b*ampV(J)+(ampV(I)*b*cos(angV(I)-angV(J)))/k;
                    H2(n,I)=-(ampV(I)*ampV(J)*b*sin(angV(I)-angV(J)))/k;
                    H2(n,J)=(ampV(I)*ampV(J)*b*sin(angV(I)-angV(J)))/k;
                else
                    H1(n,I)=ampV(J)*(g*sin(angV(I)-angV(J))+b*cos(angV(I)-angV(J)));
                    H1(n,J)=-2*ampV(J)*(b+yc)+ampV(I)*(g*sin(angV(I)-angV(J))+b*cos(angV(I)-angV(J)));
                    H2(n,I)=ampV(I)*ampV(J)*(g*cos(angV(I)-angV(J))-b*sin(angV(I)-angV(J)));
                    H2(n,J)=-ampV(I)*ampV(J)*(g*cos(angV(I)-angV(J))-b*sin(angV(I)-angV(J)));
                end
        end
    end
    if nodeRe==1                         % 删去参考节点
        H2 = H2(:,2:nbus);
    else
        H2 = [H2(:,1:(nodeRe-1)) H2(:,(nodeRe+1):nbus)];
    end
    H = [H1 H2];
end
