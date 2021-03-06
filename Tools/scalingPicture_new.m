clear
fileName_in = sprintf('results/confOriginal.dat');

%--------count number of sites-------------
linesInBlock = 0;
prevLineIsEmpty = false;
fp = fopen(fileName_in);

emptyLine = sprintf('\n');

string = fgets(fp);

blockNum = 1;
lineNum=1;

while ischar(string)
    
    if strcmp(emptyLine, string)
   
        if(~prevLineIsEmpty)
            N{blockNum}=lineNum - 1;
            blockNum = blockNum + 1;
            lineNum = 1;
        end
        prevLineIsEmpty = true;
    else
        lineNum = lineNum + 1;
        prevLineIsEmpty = false;
    end
    string = fgets(fp);
end
fclose(fp);
blockNum = blockNum-1;
%-----------------end----------------------

% read file
fp = fopen(fileName_in);

for i=1:blockNum
    N{i};
    B{i}=fscanf(fp, '%g %g %g', [3 N{i}]);
    B{i}=B{i}';
end
fclose(fp);

%axis ([-120 60 -160 0 -50 450]);
axis equal;
axis on;


k=2;
if k==1
fprintf('k=%d\n',k);
fprintf('N = %d\n',N{k});
   for i=2:N{k}-1
    %axis ([-xd xd -xd xd -zd zd]);
    %axis ([-120 60 -160 0 -50 450]);
    axis on;
    axis equal; 
    view([45, -45, 45]);
    %fprintf('Bi = %g\n',B{k}(i+1,3));
    h1=line([B{k}(i,1) B{k}(i+1,1)],[B{k}(i,2) B{k}(i+1,2)],[B{k}(i,3) B{k}(i+1,3)],'LineStyle','-','LineWidth',3,'Color',[.0 .45 .7]); 
    h2=stem3(B{k}(i,1),B{k}(i,2),B{k}(i,3),'LineStyle','none','Marker','o', 'MarkerFaceColor',[.9 0 0],'MarkerEdgeColor', [.3 0 0],'MarkerSize',10);
    hold on;
    end
%set(h2,'Visible','off')
axis equal; 
line([B{k}(1,1) B{k}(2,1)],[B{k}(1,2) B{k}(2,2)],[B{k}(1,3) B{k}(2,3)],'LineStyle','-','LineWidth',3,'Color',[.0 .4 .6]);
line([B{k}(2,1) B{k}(3,1)],[B{k}(2,2) B{k}(3,2)],[B{k}(2,3) B{k}(3,3)],'LineStyle','-','LineWidth',3,'Color',[.0 .4 .6]);
end

if k==2
fprintf('k=%d\n',k);
fprintf('%d\n',N{k}-1);
    for i=2:N{k}-1    
    %axis ([-xd xd -xd xd -zd zd]);
    %axis ([-120 60 -160 0 -50 450]);
    axis on;
    axis equal; 
    view([45, -45, 45]);
    h1=line([B{k}(i,1) B{k}(i+1,1)],[B{k}(i,2) B{k}(i+1,2)],[B{k}(i,3) B{k}(i+1,3)],'LineStyle','-','LineWidth',3,'Color',[.0 .0 .0]); 
    h2=stem3(B{k}(i,1),B{k}(i,2),B{k}(i,3),'LineStyle','none','Marker','o', 'MarkerFaceColor',[.0 .9 0],'MarkerEdgeColor', [.0 .3 0],'MarkerSize',10);
    hold on;
    end
%set(h2,'Visible','off')   
line([B{k}(1,1) B{k}(2,1)],[B{k}(1,2) B{k}(2,2)],[B{k}(1,3) B{k}(2,3)],'LineStyle','-','LineWidth',3,'Color',[.0 .0 .0]);
line([B{k}(2,1) B{k}(3,1)],[B{k}(2,2) B{k}(3,2)],[B{k}(2,3) B{k}(3,3)],'LineStyle','-','LineWidth',3,'Color',[.0 .0 .0]);
end

if k==6
fprintf('k=%d\n',k);
fprintf('%d\n',N{k}-1);
    for i=2:N{k}-1    
    %axis ([-xd xd -xd xd -zd zd]);
    %axis ([-120 60 -160 0 -50 450]);
    axis on;
    axis equal; 
    view([45, -45, 45]);
    h1=line([B{k}(i,1) B{k}(i+1,1)],[B{k}(i,2) B{k}(i+1,2)],[B{k}(i,3) B{k}(i+1,3)],'LineStyle','-','LineWidth',3,'Color',[.9 .7 .0]); 
    h2=stem3(B{k}(i,1),B{k}(i,2),B{k}(i,3),'LineStyle','none','Marker','o', 'MarkerFaceColor',[.0 .1 .9],'MarkerEdgeColor', [.0 .1 .3],'MarkerSize',10);
    hold on;
    end
%set(h2,'Visible','off')   
line([B{k}(1,1) B{k}(2,1)],[B{k}(1,2) B{k}(2,2)],[B{k}(1,3) B{k}(2,3)],'LineStyle','-','LineWidth',3,'Color',[.9 .7 .0]);
line([B{k}(2,1) B{k}(3,1)],[B{k}(2,2) B{k}(3,2)],[B{k}(2,3) B{k}(3,3)],'LineStyle','-','LineWidth',3,'Color',[.9 .7 .0]);
end
