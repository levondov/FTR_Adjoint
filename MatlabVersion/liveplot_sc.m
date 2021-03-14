hh=figure('units','pixels','position',[100,100,750,800]); 
%axis tight manual % this ensures that getframe() returns a consistent size
filename = 'ex_diff.gif';
i=1;
sp = 3;
saveGif = 0;
cc=1;
ccc=1;

Nt = length(yy);

%crts = [0,linspace(100e-6,0.75e-3,49)];
if length(z) > length(yy{1});
   z=z(1:100:end); 
end

figure(hh);
while i <= Nt
    
%% Now plot everything to compare

    subplot(2,2,1); cla(); hold on;
    plot(z,yy{i}(:,1),'k','linewidth',3);
    plot(z,yy{i}(:,2),'b','linewidth',3);
    plot(z,yy{i}(:,3),'r','linewidth',3);
    legend('q+','q-','qx','location','best');
    %ylim([-5e-5,5e-5]);
    ylabel('Q+,Q-,Qx [m^2]'); xlabel('z [m]'); grid on;

    subplot(2,2,2); cla(); hold on;
    plot(z,yy{i}(:,4),'k','linewidth',3);
    plot(z,yy{i}(:,5),'b','linewidth',3);
    plot(z,yy{i}(:,6),'r','linewidth',3);
    legend('p+','p-','px','location','best');
    %ylim([-5e-5,5e-5]);
    ylabel('P+,P-,Px [m]'); xlabel('z [m]'); grid on;

    subplot(2,2,3); cla(); hold on;
    plot(z,yy{i}(:,7),'k','linewidth',3);
    plot(z,yy{i}(:,8),'b','linewidth',3);
    plot(z,yy{i}(:,9),'r','linewidth',3);
    legend('e+','e-','ex','location','best');
    %ylim([-2e-5,3.5e-5]);
    ylabel('E+,E-,Ex '); xlabel('z [m]'); grid on;

    subplot(2,2,4); cla(); hold on;
    plot(z,yy{i}(:,10),'k','linewidth',3);
    %plot(z,y(:,11),'r','linewidth',3);
    legend('L','location','best');
    %ylim([-3e-5,5e-6]);
    ylabel('L'); xlabel('z [m]'); grid on;
    
    tstring = [num2str(crts(i)*1e3),' mA'];
    subplot(2,2,1); title(tstring);
    subplot(2,2,2); title(tstring);
  
    
    if ~saveGif
        k = waitforbuttonpress;
    else
        if i==1 && cc==1
            gif(filename,'DelayTime',0.2,'frame',hh);
        else
            gif;
        end
    end
    % 28 leftarrow
    % 29 rightarrow
    % 30 uparrow
    % 31 downarrow
    value = double(get(gcf,'CurrentCharacter'));
    if value == 28
        i = i-sp;
    elseif value == 29
        i = i+sp;
    elseif value == 30
        i = 1;
    elseif value == 31
        i = Nt;
    else
        if cc < 10
            cc=cc+1;
        else
            i=i+sp;
        end
    end
    
    if i < 1
        i = 1;
    end
    if i > Nt
        if ~saveGif
            i = Nt;
        else
            if ccc < 10
                ccc = ccc+1;
                i=Nt;
            end
        end
    end
end


