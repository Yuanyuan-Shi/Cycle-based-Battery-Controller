clc
clear
rng(3);
r = rand(10,1);
r(r>.5)=1;
r(r<=.5)=-1;
ss = -cumsum(r);
ss = [0;ss(1:(end-1))];
ext=sig2ext(ss);
a=rainflow(ext,1);
[m n]=size(a);

% if n>100,
%     button = questdlg(['Rainflow found ' num2str(sum(a(3,:))) ' cycles! Do you want to continue?'],...
%         'Continue Operation','Yes','No','No');
%     if strcmp(button,'No')
%         error('Function aborted by user.')        
%     end
% end

col='ymcrgb';
plot(0:length(ext)-1,ext,'k.:')
hold on
wyk=0:0.05:1;
for c=1:n,
    colnr=rem(c-1,6)+1;
    
    nr1=round(a(4,c)+1);
    nr2=round(a(4,c)+1+a(5,c)*a(3,c));
    if a(3,c)==1.0,
        if ext(nr1)<ext(nr1+1),
            plot(wyk.*a(5,c)+a(4,c),cos(pi+wyk.*2*pi)*a(1,c)+a(2,c),col(colnr))
            text(a(4,c),a(2,c)-a(1,c),[int2str(c) '. Cycle, up'],...
                'Color',col(colnr),'VerticalAlignment','top')
        else
            plot(wyk.*a(5,c)+a(4,c),cos(   wyk.*2*pi)*a(1,c)+a(2,c),col(colnr))
            text(a(4,c),a(2,c)+a(1,c),[int2str(c) '. Cycle, down'],...
                'Color',col(colnr),'VerticalAlignment','bottom')
        end
    else
        if ext(nr1)>ext(nr2),
            plot(wyk.*a(5,c)*0.5+a(4,c),cos(   wyk.*pi)*a(1,c)+a(2,c),col(colnr))
            text(a(4,c),a(2,c)+a(1,c),[int2str(c) '. Half-cycle, down'],...
                'Color',col(colnr),'VerticalAlignment','bottom')
        else
            plot(wyk.*a(5,c)*0.5+a(4,c),cos(pi+wyk.*pi)*a(1,c)+a(2,c),col(colnr))
            text(a(4,c),a(2,c)-a(1,c),[int2str(c) '. Half-cycle, up'],...
                'Color',col(colnr),'VerticalAlignment','top')
        end
    end
end
xlabel('peaks, counted from 0')
ylabel('value')
title('Rainflow cycles extracted from signal')
legend('peaks from signal',0)
hold off

disp('Row 1: amplitude')
disp('Row 2: mean')
disp('Row 3: number of cycles (cycle or half cycle)')
disp('Row 4: begin time of extracted cycle or half cycle')
disp('Row 5: period of a cycle')   
disp(a)