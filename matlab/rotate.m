function [uppermodrot,lowermodrot]=rotate(alfa)
    load airfoilrescalled.mat
    xo=0.25;
    yo=0;
    [ma,na]=size(uppermod);
    [mb,nb]=size(lowermod);
    alfa=-alfa*(2*pi/360);
    for i=1:ma
        uppermodrot(i,1)=xo+cos(alfa)*(uppermod(i,1)-xo)-sin(alfa)*(uppermod(i,2)-yo);
        uppermodrot(i,2)=yo+sin(alfa)*(uppermod(i,1)-xo)+cos(alfa)*(uppermod(i,2)-yo);
    end
    for i=1:mb
        lowermodrot(i,1)=xo+cos(alfa)*(lowermod(i,1)-xo)-sin(alfa)*(lowermod(i,2)-yo);
        lowermodrot(i,2)=yo+sin(alfa)*(lowermod(i,1)-xo)+cos(alfa)*(lowermod(i,2)-yo);
    end
    alfa=-alfa*(360/2/pi);
    str=['airfoil_rescalled_rotated' num2str(alfa) '.dat'];
    fid=fopen(str, 'w');
    fprintf(fid,'%d \n', ma+mb);
    for i=1:ma
        fprintf(fid,'%f %f %d \n', uppermodrot(i,1), uppermodrot(i,2), 0);
    end
    for i=1:mb
        fprintf(fid,'%f %f %d \n', lowermodrot(i,1), lowermodrot(i,2), 0);
    end
    plot(lowermod(:,1),lowermod(:,2));
    hold on
    plot(uppermod(:,1),uppermod(:,2));
    plot(lowermodrot(:,1),lowermodrot(:,2),'*');
    plot(uppermodrot(:,1),uppermodrot(:,2),'o');
    fclose(fid);
end