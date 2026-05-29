clear all; close all;
%folder = ''; % Give the path to the file.  
 
fileID='Example_data_over_time_WT_condition.dat'; % Main data-file name

hex = dlmread(fileID);  

 T=hex(:,1); % 1st column contains all the Time counter
 cell = hex(:,2); % 2nd column contains cell index
 Vx=hex(:,4); % 4th column is bead pos (X)
 Vy=hex(:,5); % 5th column is bead pos (Y)
 a = 24.0; % Ellipse major axis length
 b = 12.0; % Ellipse minor axis length
 
ymax=max(Vy(:));
%ymax=box;
ymin=min(Vy(:));
%ymin=0.0;
%xmax=max(Vx(:));
xmax=a;
xmin=min(Vx(:));
%xmin=0.0;

cellStart = min(cell); 
cellEnd = max(cell); 

tmin = min(T); % Time(iteration) at beginning.
tmax= 500000;  % Time(iteration) upto which we want the movie

dt = 4000; % How often would you like to capture the frame.
dcell=1;

 %% Define movie name etc. 
 
 figure('units','pixels','position',[0 0 1920 1080])
% figure(1), % Shall open an empty figure window. Adjust it to the size you want. 
 % Ideally, you should have a fixed figure with specified position and size, but 
 % I find that very constrained. So, I just open an empty figure window, and adjust 
 % it according to my need. Don't touch it until the code completes. It'll interrupt 
 % the program, and the movie won't be produced. I sometime use this method to kill a 
 % process and redo the movie to my like.   
 
% clear movObj

movObj = VideoWriter('polyfill_WT.avi','Uncompressed AVI'); % Give whatever name you want to give. Can be automated. 
movObj.FrameRate = 10; % Change to whatever value you want. 
Quality = 100; % Any number between 0-100. 
open(movObj) 

 
 %% The best part: Frame capture 
ymax_oocyte = b*sqrt(0.25-((0.75*a-a/2.0)/a)^2);
y_oocyte = linspace(ymax_oocyte,-ymax_oocyte,50); %% oocyte left boundary
x_oocyte = zeros(1,50) + 0.75*a ;
theta=0 ;
x_oocyte_upper = [];
x_oocyte_lower = [];
y_oocyte_upper = [];
y_oocyte_lower = [];


 for j=1:500
         x1(j)=(a/2.0)*(1+cos(theta));
         y1(j)=(b/2.0)*sin(theta);
         theta = theta + 2*pi/500;
         
         if x1(j)>= 0.75*a && y1(j)>=0
             x_oocyte_upper = [x_oocyte_upper, x1(j)];
             y_oocyte_upper = [y_oocyte_upper, y1(j)];
         elseif x1(j)> 0.75*a && y1(j)<0
             x_oocyte_lower = [x_oocyte_lower, x1(j)];
             y_oocyte_lower = [y_oocyte_lower, y1(j)];
         end    
 end 
 
 x_oocyte = [x_oocyte, x_oocyte_lower, x_oocyte_upper]; 
 y_oocyte = [y_oocyte, y_oocyte_lower, y_oocyte_upper];
 
 for t = tmin:dt:tmax 

     plot(x1,y1,'Color',[165, 42, 42]./255,'LineWidth',12); 
     hold on;
     patch(x_oocyte,y_oocyte,[0.8 0.8 0.8],'EdgeColor',[165, 42, 42]./255,'LineStyle',"--",'LineWidth',2);
     plot(x_oocyte_lower,y_oocyte_lower,'black','LineWidth',12);
     plot(x_oocyte_upper,y_oocyte_upper,'black','LineWidth',12);
     
     
 for cellindex = cellStart:dcell:cellEnd 
     
	id = cell == cellindex & T==t; % All data index for time t and for cell index 'i'. 
    
	x = Vx(id);
    y = Vy(id);
    
    if cellindex>=13 % Polar cells
    patch(x,y,'magenta','EdgeColor','black','LineWidth',1.0)  % patch function fills the polygons; 
    hold on;
    elseif cellindex<=6 % Nurse cells
    patch(x,y,[0.8 0.8 0.8],'EdgeColor','black','LineWidth',7) % patch function fills the polygons; grey color code; 
    hold on;
    elseif cellindex==7 % Border cell
    patch(x,y,[0,109,44]./255)  % patch function fills the polygons; 
    hold on;         
    plot(x,y,'black-','LineWidth',1);
    
    elseif cellindex==8 % Border cell
    patch(x,y,[161,217,155]./255)  % patch function fills the polygons; grey color code;    
    hold on;          
    plot(x,y,'black-','LineWidth',1);
    
    elseif cellindex==9 % Border cell
    patch(x,y,[49,163,84]./255)  % patch function fills the polygons; grey color code;    
    hold on;        
    plot(x,y,'black-','LineWidth',1);
    
    elseif cellindex==10 % Border cell
    patch(x,y,[199,233,192]./255)  % patch function fills the polygons; grey color code;    
    hold on;        
    plot(x,y,'black-','LineWidth',1);
    
    elseif cellindex==11 % Border cell
    patch(x,y,[116,196,118]./255)  % patch function fills the polygons; grey color code;    
    hold on;         
    plot(x,y,'black-','LineWidth',1);
     
    elseif cellindex==12 % Border cell
    patch(x,y,[237,248,233]./255)     
    hold on;         
    plot(x,y,'black-','LineWidth',1);
    end    
    
    title('Simulation Video: WT','Fontsize', 20)
    axis equal
    axis off;
	xlim([0 24]);
    ylim([-6.3 6.3]);

 end

    hold off;    
    M(t)=getframe;  
	frame = getframe(gcf); % Capture the frame
    writeVideo(movObj,frame); % Add the frame to the movie.
 end

 close(movObj); % Close the movie variable and write the file. 
