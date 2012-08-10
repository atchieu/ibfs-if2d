function [uppermod,lowermod]=transform(upper,lower,d,inv)
    xal=lower(:,1);
    xau=upper(:,1);
    yal=lower(:,2);
    yau=upper(:,2);
    i=1;
    xmodu(1)=0;
    xmodl(1)=0;
    check1=0;
    check2=0;
    bdynou(1)=0;
    bdynol(1)=0;
    while check1==0 || check2==0
        if check1==0
            pom1=xmodu(i)+abs(fsolve( @(dx)nextnode(dx,xmodu(i),d,xau,yau),xmodu(i)));
            if pom1<1
                xmodu(i+1)=pom1;
                bdynou(end+1)=0;
            else
                check1=1;
            end
        end
        if check2==0    
            pom2=xmodl(i)+abs(fsolve( @(dx)nextnode(dx,xmodl(i),d,xal,yal),xmodl(i)));
            if pom2<1
                xmodl(i+1)=pom2;
                bdynol(end+1)=0;
            else
            check2=1;
            end
        end
        i=i+1;
    end
    %xmodl(1)=[];
    %bdynol(1)=[]
    ymodu=airfoil(xmodu, xau,yau);
    ymodl=airfoil(xmodl, xal,yal);
    uppermod=[xmodu' ymodu'];
    lowermod=[xmodl' ymodl'];
    [m,n]=size(uppermod);
    j=1;
    if inv==1
        for i=m:-1:1
            uppermodinv(j,:)=uppermod(i,:);
            j=j+1;
        end
        uppermod=uppermodinv;
    end
    uppermod=[xmodu' ymodu' bdynou'];
    lowermod=[xmodl' ymodl' bdynol'];
    save('airfoilrescalled.mat','uppermod', 'lowermod')
    fid=fopen('airfoilrescalled.dat', 'w');
    if fid
        f=1;
    end
    plot(uppermod(:,1),uppermod(:,2),'o');
    hold on
    plot(lowermod(:,1),lowermod(:,2),'*');
    [ma,na]=size(uppermod);
    [mb,nb]=size(lowermod);
    fprintf(fid,'%d \n', ma+mb);
    for i=1:ma
        fprintf(fid,'%f %f %d \n', uppermod(i,1), uppermod(i,2), uppermod(i,3));
    end
    for i=1:mb
        fprintf(fid,'%f %f %d \n', lowermod(i,1), lowermod(i,2), lowermod(i,3));
    end
    fclose(fid);
    

end