clear all;
clc;
I=2;
KJ = 10^3;

u=[-200,-450;100,-280;0,200;-350,-150;350,130;-100,280; ]; 
zz=[280,140;120,400;240,-80;-50,-250;-260,-30;-50,0];
c=[300,20;-80,200;-200,0];

w=[7; 13; 4; 16; 10; 8;]/KJ;%30; 20; 35
v=[10; 10; 10; 10; 10; 10]/KJ;%20; 10; 15 

final_rate=zeros(100,1);
energy=zeros(100,1);

load("u_data.mat")
load("zz_data.mat")

initial=[30;30]/KJ;

mu=100;
H=50;

for iter=1:100

    Pi=1; N=100; K=6; M=6; Vmax=20; Vz=10; Tp=5; Td=5; R=[100,70,50,50,50,40]; e=1e-1; dmin=10;
    Hmin1=20; Hmax1=20;
    Hmin2=25; Hmax2=25; 
    Wmax=[60/KJ,60/KJ;];
    UAV=2/KJ; 
    P0=79.85/KJ; Ut=120; d=0.6; rho=1.225/KJ; s0=0.05; A0=0.503; v0=4.03; g=9.8; P1=1.1/sqrt(2*rho*A0); %energy parameter
    
    qstr=[-390 -410;-410 -390]; qdst=[410 390;390 410];
    zstr=[Hmin1;Hmin2]; zdst=[Hmax1;Hmax2];
    
    b=zeros(N,K,I); a=zeros(N,K,I);  t_pre=0.1.*ones(N,I); 
    Emax=[1500,1500];   
    
    del_=ones(N,I);
    x0=linspace(qstr(1,1),qdst(1,1),N)';
    y0=linspace(qstr(1,2),qdst(1,2),N)';
    
    x1=linspace(qstr(2,1),qdst(2,1),N)';
    y1=linspace(qstr(2,2),qdst(2,2),N)';
    q1=[x0,y0];
    q2=[x1,y1];
    q(:,:,1)=q1;
    q(:,:,2)=q2;
    
    z0=linspace(zstr(1),zdst(1),N)';
    z1=linspace(zstr(2),zdst(2),N)';
    z=[z0,z1];
        
    for i=1:I
        for n=2:N
            qnorm(n,i)=norm(q(n,:,i)-q(n-1,:,i));
        end
    end
    
    
    ZZ=3;
    NFZ=zeros(N,ZZ,2);
    for m=1:ZZ
        for n=1:N
            NFZ(n,m,:)=[c(m,1)+R(m)*cos(2*pi*(n-1)/(N-1)),c(m,2)+R(m)*sin(2*pi*(n-1)/(N-1))];
        end
    end
    
    %     load("user6.mat")
    %     mu=100;
    %     u=u_data(:,:,iter);
    %     zz=zz_data(:,:,iter);

    eta1=zeros(100,1);
    eta2=zeros(100,1);
    count=1;

    tic
    while(1)
    
        for i=1:I
            for n=1:N
                Wr(n,i)=sum(sum(a(1:n,:,i)*w)-sum(b(1:n,:,i)*v))+((initial(i)));
            end
        end
        
        for i=1:I
            for n=2:N
                qnorm(n,i)=norm(q(n,:,i)-q(n-1,:,i));
            end
        end
    
        cvx_begin quiet
        variables A(N,K,I) AA(N,K,I) B(N,M,I) BB(N,M,I) Q(N,2,I) del(N,I) Z(N,I) t(N,I) tt(N,I)
        expressions W(N,I) 
    
        subject to
           
        0<=t
        0<=del
        0<=Z
    
        B.*(1-2*b)+b.^2<=BB
        0<=BB
        0<=B<=1
        
        A.*(1-2*a)+a.^2<=AA
        0<=AA
        0<=A<=1

        Q(1,:,1)==qstr(1,:); Q(N,:,1)==qdst(1,:); 
        Q(1,:,2)==qstr(2,:); Q(N,:,2)==qdst(2,:); 
        Z(1,:)==zstr'; Z(N,:)==zdst'; 

        for i=1:I
            norms([Z(2:N,i)-Z(1:N-1,i) Q(2:N,1,i)-Q(1:N-1,1,i) Q(2:N,2,i)-Q(1:N-1,2,i)],2,2)<=Vmax*del(2:N,i)
        end

        for i=1:I
            norms(Z(2:N,i)-Z(1:N-1,i),2,2)<=Vz*del(2:N,i)
        end
    
        for i=1:I
            for m=1:ZZ
                -((q(:,1,i)-c(m,1)).^2+(q(:,2,i)-c(m,2)).^2)-2*((q(:,1,i)-c(m,1)).*(Q(:,1,i)-q(:,1,i))+(q(:,2,i)-c(m,2)).*(Q(:,2,i)-q(:,2,i)))<=-(R(m)^2+Vmax.^2*del(:,i).^2/4)
                -((q(1:N-1,1,i)-c(m,1)).^2+(q(1:N-1,2,i)-c(m,2)).^2)-2*((q(1:N-1,1,i)-c(m,1)).*(Q(1:N-1,1,i)-q(1:N-1,1,i))+(q(1:N-1,2,i)-c(m,2)).*(Q(1:N-1,2,i)-q(1:N-1,2,i)))<=-(R(m)^2+Vmax.^2*del(2:N,i).^2/4)
            end
        end

        (1-sum(B(:,:,1),2)-sum(A(:,:,1),2))*Hmin1<=Z(:,1)<=Hmax1
        (1-sum(B(:,:,2),2)-sum(A(:,:,2),2))*Hmin2<=Z(:,2)<=Hmax2

        for i=1:I
            for k=1:K
                A(:,k,i)<=-(H^2+Pi^2)./(z(:,i).^2+(q(:,1,i)-u(k,1)).^2+(q(:,2,i)-u(k,2)).^2+H^2).^2.*(((Q(:,1,i)-u(k,1)).^2+(Q(:,2,i)-u(k,2)).^2+Z(:,i).^2)-((q(:,1,i)-u(k,1)).^2+(q(:,2,i)-u(k,2)).^2+z(:,i).^2))+(H^2+Pi^2)./(z(:,i).^2+(q(:,1,i)-u(k,1)).^2+(q(:,2,i)-u(k,2)).^2+H^2)
            end
        end
    
        for i=1:I
            for k=1:K
                A(2:N,k,i)<=-(H^2+Pi^2)./((q(1:N-1,1,i)-u(k,1)).^2+(q(1:N-1,2,i)-u(k,2)).^2+H^2).^2.*(((Q(1:N-1,1,i)-u(k,1)).^2+(Q(1:N-1,2,i)-u(k,2)).^2)-((q(1:N-1,1,i)-u(k,1)).^2+(q(1:N-1,2,i)-u(k,2)).^2))+(H^2+Pi^2)./((q(1:N-1,1,i)-u(k,1)).^2+(q(1:N-1,2,i)-u(k,2)).^2+H^2)
            end
        end
    
        for i=1:I
            for k=1:K
                A(1:N-1,k,i)<=-(H^2+Pi^2)./((q(2:N,1,i)-u(k,1)).^2+(q(2:N,2,i)-u(k,2)).^2+H^2).^2.*(((Q(2:N,1,i)-u(k,1)).^2+(Q(2:N,2,i)-u(k,2)).^2)-((q(2:N,1,i)-u(k,1)).^2+(q(2:N,2,i)-u(k,2)).^2))+(H^2+Pi^2)./((q(2:N,1,i)-u(k,1)).^2+(q(2:N,2,i)-u(k,2)).^2+H^2)
            end
        end
    
        sum(sum(A,3))<=1   
        sum(sum(sum(A)))>=K
    
        sum(A,2)<=1
        sum(B,2)<=1
        
        for i=1:I
            for m=1:M
                B(:,m,i)<=-(H^2+Pi^2)./(z(:,i).^2+(q(:,1,i)-zz(m,1)).^2+(q(:,2,i)-zz(m,2)).^2+H^2).^2.*(((Q(:,1,i)-zz(m,1)).^2+(Q(:,2,i)-zz(m,2)).^2+Z(:,i).^2)-((q(:,1,i)-zz(m,1)).^2+(q(:,2,i)-zz(m,2)).^2+z(:,i).^2))+(H^2+Pi^2)./(z(:,i).^2+(q(:,1,i)-zz(m,1)).^2+(q(:,2,i)-zz(m,2)).^2+H^2)
            end
        end
    
        for i=1:I
            for m=1:M
                B(2:N,m,i)<=-(H^2+Pi^2)./((q(1:N-1,1,i)-zz(m,1)).^2+(q(1:N-1,2,i)-zz(m,2)).^2+H^2).^2.*(((Q(1:N-1,1,i)-zz(m,1)).^2+(Q(1:N-1,2,i)-zz(m,2)).^2)-((q(1:N-1,1,i)-zz(m,1)).^2+(q(1:N-1,2,i)-zz(m,2)).^2))+(H^2+Pi^2)./((q(1:N-1,1,i)-zz(m,1)).^2+(q(1:N-1,2,i)-zz(m,2)).^2+H^2)
            end
        end
    
        for i=1:I
            for m=1:M
                B(1:N-1,m,i)<=-(H^2+Pi^2)./((q(2:N,1,i)-zz(m,1)).^2+(q(2:N,2,i)-zz(m,2)).^2+H^2).^2.*(((Q(2:N,1,i)-zz(m,1)).^2+(Q(2:N,2,i)-zz(m,2)).^2)-((q(2:N,1,i)-zz(m,1)).^2+(q(2:N,2,i)-zz(m,2)).^2))+(H^2+Pi^2)./((q(2:N,1,i)-zz(m,1)).^2+(q(2:N,2,i)-zz(m,2)).^2+H^2)
            end
        end
    
        sum(sum(B,3))<=1
        10/KJ*(sum(sum(B)))>=reshape(initial,[1,1,2])
        
        for i=1:I
            for n=1:N
                W(n,i)=sum(sum(A(1:n,:,i)*w)-sum(B(1:n,:,i)*v))+((initial(i)));
            end
        end
           
        for i=1:I
            W(:,i)<=Wmax(:,i)
        end
    
        for i=1:I
            idx1 = 2:N;
            idx2 = 1:N-1;
            pow_pos( quad_over_lin( del(idx1,i) ,t(idx1,i),0),2 ) <= t_pre(idx1,i).^2 + 2*t_pre(idx1,i).*(t(idx1,i)-t_pre(idx1,i)) - ( (q(idx1,1,i)-q(idx2,1,i)).^2+(q(idx1,2,i)-q(idx2,2,i)).^2 )/v0^2 + 2*( (q(idx1,1,i)-q(idx2,1,i)).*(Q(idx1,1,i)-Q(idx2,1,i)) + (q(idx1,2,i)-q(idx2,2,i)).*(Q(idx1,2,i)-Q(idx2,2,i)) )/v0^2  % 45 -> (47)
            geo_mean([tt(idx1,i),del(idx1,i),del(idx1,i)],2) >= norms(Q(idx1,:,i)-Q(idx2,:,i),2,2)
    
        end
        
        for i=1:I
            idx1 = 2:N;
            idx2 = 1:N-1;
        
            sum(P0*( del(idx1,i)+3/Ut^2*(quad_over_lin( qnorm(idx1,i) ,del(idx1,i),0)) ) +...
                1/2*d*rho*s0*A0*tt(idx1,i) + P1.*((  square_pos( ((W(idx1,i)+UAV)*g).^(3/2)+t(idx1,i) ) - ( ( (Wr(idx1,i)+UAV)*g).^3 + 3*((Wr(idx1,i)+UAV)*g).^2.*( ((W(idx1,i)+UAV)*g) - ((Wr(idx1,i)+UAV)*g) ) ) - ( (t_pre(idx1,i)).^2+2*t_pre(idx1,i).*(t(idx1,i)-t_pre(idx1,i)) )  )./2) + ...
                ((  square_pos( ((W(idx1,i)+UAV)*g)+abs(Z(idx1,i)-Z(idx2,i)) ) - ( ( (Wr(idx1,i)+UAV)*g).^2 + 2*((Wr(idx1,i)+UAV)*g).*( ((W(idx1,i)+UAV)*g) - ((Wr(idx1,i)+UAV)*g) ) ) - ( -(z(idx1,i)-z(idx2,i)).^2 + 2*(z(idx1,i)-z(idx2,i)).*( (Z(idx1,i)-Z(idx2,i)) ) ) )./2))<=Emax(i)
        
        end
    
        minimize sum(sum(del))+mu*(sum(sum(sum(AA)))+sum(sum(sum(BB))))
    
        cvx_end
        eta1(count)=sum(sum(del))+mu*(sum(sum(sum(AA)))+sum(sum(sum(BB))));
        eta2(count)=sum(sum(del));

        if strcmp(cvx_status,'Infeasible')~=0
            cvx_status
            break
        end

        if count>=100
            break
        end
    
        count=count+1;

        mu=min(mu*1.5,100000);

        b=B; a=A; q=Q; z=Z; t_pre=t;
    
        for i=1:I
            for n=2:N
                qnorm(n,i)=norm(q(n,:,i)-q(n-1,:,i));
                znorm(n,i)=norm(z(n,i)-z(n-1,i));
            end
        end
    
        for i=1:I
            idx1 = 2:N;
            idx2 = 1:N-1;
            Energy(i) = sum(P0*( del(idx1,i)+3/Ut^2*(quad_over_lin( qnorm(idx1,i) ,del(idx1,i),0)) ) +...
                1/2*d*rho*s0*A0*tt(idx1,i) + P1.*((  square_pos( ((W(idx1,i)+UAV)*g).^(3/2)+t(idx1,i) ) - ( ( (Wr(idx1,i)+UAV)*g).^3 + 3*((Wr(idx1,i)+UAV)*g).^2.*( ((W(idx1,i)+UAV)*g) - ((Wr(idx1,i)+UAV)*g) ) ) - ( (t_pre(idx1,i)).^2+2*t_pre(idx1,i).*(t(idx1,i)-t_pre(idx1,i)) )  )./2) + ...
                ((  square_pos( ((W(idx1,i)+UAV)*g)+abs(Z(idx1,i)-Z(idx2,i)) ) - ( ( (Wr(idx1,i)+UAV)*g).^2 + 2*((Wr(idx1,i)+UAV)*g).*( ((W(idx1,i)+UAV)*g) - ((Wr(idx1,i)+UAV)*g) ) ) - ( -(z(idx1,i)-z(idx2,i)).^2 + 2*(z(idx1,i)-z(idx2,i)).*( (Z(idx1,i)-Z(idx2,i)) ) ) )./2));
        end
        
        for i=1:I
            idx1 = 2:N;
            idx2 = 1:N-1;
            Energy_real(i) = sum(P0*( del(idx1,i)+3/Ut^2*qnorm(idx1,i).^2./del(idx1,i) ) +...
            1/2*d*rho*s0*A0*(qnorm(idx1,i)).^3./del(idx1,i).^2 + P1.*( ((W(idx1,i)+UAV)*g).^(3/2).*(sqrt(del(idx1,i).^4+qnorm(idx1,i).^4./(4*v0^4))-qnorm(idx1,i).^2./(2*v0^2)).^(1/2) ) + ...
            ((W(idx1,i)+UAV)*g).*znorm(idx1,i));
        end
        
        if(abs(sum(sum(del_))-sum(sum(del)))+sum(sum(sum(AA)))+sum(sum(sum(BB)))<e)
            break;
        end
        del_=del;
    
        drawnow
        clf
        hold on
        xlim([-500 500])
        ylim([-500 500])
        scatter(qstr(:,1),qstr(:,2),'*')
        scatter(qdst(:,1),qdst(:,2),'*')
        scatter(zz(1:K,1),zz(1:K,2))
        scatter(u(1:K,1),u(1:K,2))
    
        
        for i=1:I
            plot3(Q(:,1,i),Q(:,2,i),Z(:,i))
        end
    
        for m=1:ZZ
            plot(NFZ(:,m,1),NFZ(:,m,2));
        end
    
        for n=1:N
            for k=1:K
                for i=1:I
                    if b(n,k,i)<0
                        b(n,k,i)=0;
                    end
                    if a(n,k,i)<0
                        a(n,k,i)=0;
                    end
                end
            end
        end
    end

    toc

    if strcmp(cvx_status,'Infeasible')~=0 || count>=100
        cvx_status
    else
        final_rate(iter,1)=sum(sum(del));
    end

    save("scenario2_user6","final_rate")

end

