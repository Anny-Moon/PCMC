%  Copyright 2017 Anna Sinelnikova
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%       http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

function [] = makeThermalizationMovie(increment)

% default argument for increment
if nargin < 1
    increment =1 ;
end

fileNameIn = sprintf('results/Configurations/0confR.dat');
fileNameOut = sprintf('thermalization.avi');

%--------count number of blocks and sites-------------
linesInBlock = 0;
prevLineIsEmpty = false;
fp = fopen(fileNameIn);

emptyLine = sprintf('\n');

string = fgets(fp);

blockNum = 1;
lineNum=1;

while ischar(string)
    
    if strcmp(emptyLine, string)
   
        if(~prevLineIsEmpty)
            N{blockNum}=lineNum - 1;% fill matrix with number of lines in the blok
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
%---------------------end------------------------

% read file
fp = fopen(fileNameIn);

for i=1:blockNum
    N{i};
    B{i}=fscanf(fp, '%g %g %g', [3 N{i}]);
    B{i}=B{i}';
end
fclose(fp);


% write video
writerObj = VideoWriter(fileNameOut);% avi name
writerObj.FrameRate = 2;% frames per second
%writeObj.Quaity = 100;
%writeObj.Height=[2000,2000];

open(writerObj);
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','opengl','Position',[0 0 500 500]);

%set(gcf,'color','white')

% write the first frame (the original chain)
k=1;
fprintf('Number of monomers in initial chain: %d\n',N{k}-1);
for i=2:N{k}-1    
    axis on;
    axis equal; 
    view([60, 60, 45]);
    h1=line([B{k}(i,1) B{k}(i+1,1)],[B{k}(i,2) B{k}(i+1,2)],[B{k}(i,3) B{k}(i+1,3)],'LineStyle','-','LineWidth',3,'Color',[.0 .45 .7]); 
    h2=stem3(B{k}(i,1),B{k}(i,2),B{k}(i,3),'LineStyle','none','Marker','o', 'MarkerFaceColor',[.9 0 0],'MarkerEdgeColor', [.3 0 0],'MarkerSize',10);
    hold on;
end 

view([60, 60, 45]);
line([B{k}(1,1) B{k}(2,1)],[B{k}(1,2) B{k}(2,2)],[B{k}(1,3) B{k}(2,3)],'LineStyle','-','LineWidth',3,'Color',[.0 .4 .6]);
line([B{k}(2,1) B{k}(3,1)],[B{k}(2,2) B{k}(3,2)],[B{k}(2,3) B{k}(3,3)],'LineStyle','-','LineWidth',3,'Color',[.0 .4 .6]);
axis equal; 
hold on;

frame = getframe(gcf);
writeVideo(writerObj,frame);

% write all the others frames (chains after scaling)
for k=increment:increment:blockNum
 hold off   

    for i=2:N{k}-1    
        view([60, 60, 45]);
        h1=line([B{k}(i,1) B{k}(i+1,1)],[B{k}(i,2) B{k}(i+1,2)],[B{k}(i,3) B{k}(i+1,3)],'LineStyle','-','LineWidth',3,'Color',[.0 .45 .7]); 
        h2=stem3(B{k}(i,1),B{k}(i,2),B{k}(i,3),'LineStyle','none','Marker','o', 'MarkerFaceColor',[.9 0 0],'MarkerEdgeColor', [.3 0 0],'MarkerSize',10);
 
        axis equal;
%        axis ([xLimits yLimits zLimits]);% set axes limits as from the first frame
        axis on;
        hold on;
    end
    
    view([60, 60, 45]);
    line([B{k}(1,1) B{k}(2,1)],[B{k}(1,2) B{k}(2,2)],[B{k}(1,3) B{k}(2,3)],'LineStyle','-','LineWidth',3,'Color',[.0 .4 .6]);
    line([B{k}(2,1) B{k}(3,1)],[B{k}(2,2) B{k}(3,2)],[B{k}(2,3) B{k}(3,3)],'LineStyle','-','LineWidth',3,'Color',[.0 .4 .6]);
    
    axis equal;
%    axis ([xLimits yLimits zLimits]);% set axes limits as from the first frame
    axis on;
    fprintf('step = %d\n',k);    
    frame = getframe(gcf);
    writeVideo(writerObj,frame);

end

close(writerObj);
%close;

clear;