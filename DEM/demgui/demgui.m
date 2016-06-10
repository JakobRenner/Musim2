function demgui
f=figure('MenuBar','none','Name','Granular-Simulator','NumberTitle','off','Position',[5,50,900,600]);%,'Toolbar','figure');
defaultBackground = get(0,'defaultUicontrolBackgroundColor');
set(f,'Color',defaultBackground)

Title = uicontrol('Style','Text','String','Granular Dynamics - Simulator','Fontsize',18,'Position',[175,540,550,50],'HorizontalAlignment','center');

axes('Position',[0.89, 0.88, 0.12, 0.12]);
logo=imread('ww8-logo.png','png');
image(logo);
axis image;
axis off;

List = {'Rubber','Aluminium','Glass','Silicon','Steel','Polypropylene'};
ListText = uicontrol('Style','Text','String','Material:','Position',[10,500,100,20],'HorizontalAlignment','left','FontSize',8);  
ListBox = uicontrol('Style','popupmenu','String',List,'Position',[120,500,162,20]);

FricText = uicontrol('Style','Text','String','Friction:','Position',[10,470,100,20],'HorizontalAlignment','left','FontSize',9);
FricEdit = uicontrol('Style','Edit','String','0.3','Position',[120,470,50,20],'HorizontalAlignment','center');
FricSlider = uicontrol('Style','Slider','Position',[180,470,100,20],'CallBack', @FricSliderCallBack, 'Value',0.3,'Min',0,'Max',1);

GravText = uicontrol('Style','Text','String','Gravity:','Position',[10,440,100,20],'HorizontalAlignment','left','FontSize',9);
GravEdit = uicontrol('Style','Edit','String','9.81','Position',[120,440,50,20],'HorizontalAlignment','center');
GravSlider = uicontrol('Style','Slider','Position',[180,440,100,20],'CallBack', @GravSliderCallBack, 'Value',9.81,'Min',0,'Max',20);

FileText = uicontrol('Style','Text','String','File:','Position',[10,330,100,20],'HorizontalAlignment','left','FontSize',9);
FileEdit = uicontrol('Style', 'Edit', 'String','no file selected','Position',[10,310,150,20]);
uicontrol('Style','PushButton','String','Browse','Position',[180,310,100,20],'CallBack',@BrowseButtonPressed);


StepText = uicontrol('Style','Text','String','Time Step (s):','Position',[10,410,100,20],'HorizontalAlignment','left','FontSize',9);
StepEdit = uicontrol('Style','Edit','String','1e-6','Position',[120,410,50,20],'HorizontalAlignment','center');
StepSlider = uicontrol('Style','Slider','Position',[180,410,100,20],'CallBack', @StepSliderCallBack, 'Value',1e-6,'Min',1e-8,'Max',1e-3);

TimeText = uicontrol('Style','Text','String','Run Time (s):','Position',[10,380,100,20],'HorizontalAlignment','left','FontSize',9);
TimeEdit = uicontrol('Style','Edit','String','5','Position',[120,380,50,20],'HorizontalAlignment','center');
TimeSlider = uicontrol('Style','Slider','Position',[180,380,100,20],'CallBack', @TimeSliderCallBack, 'Value',5,'Min',0,'Max',10);

uicontrol('Style','PushButton','String','Snapshot','Position',[70,230,150,60],'CallBack',@SnapCallback);

uicontrol('Style','PushButton','String','Start','Position',[70,150,150,60],'Interruptible','on','CallBack',@StartButtonPressed);

uicontrol('Style','PushButton','String','Close','Position',[70,70,150,60],'CallBack','close');

uicontrol('Style','Text','String','Time:','Position',[50,20,80,20],'HorizontalAlignment','left','FontSize',9);
CurrentTimeText = uicontrol('Style','Text','String','0','Position',[135,20,70,20],'HorizontalAlignment','center','FontSize',9);
uicontrol('Style','Text','String','s','Position',[210,20,40,20],'HorizontalAlignment','left','FontSize',9);

ah = axes('Position',[0.35,0.075,0.6,0.8]);

    function YoungSliderCallBack(varargin)
        
        num = get(YoungSlider, 'Value');
        set(YoungEdit, 'String', num2str(num));
        
    end
    function TimeSliderCallBack(varargin)
        
        num = get(TimeSlider, 'Value');
        set(TimeEdit, 'String', num2str(num));
        
    end
    function StepSliderCallBack(varargin)
        
        num = get(StepSlider, 'Value');
        set(StepEdit, 'String', num2str(num));
        
    end
    function FricSliderCallBack(varargin)
        
        num = get(FricSlider, 'Value');
        set(FricEdit, 'String', num2str(num));
        
    end 
    function GravSliderCallBack(varargin)
        
        num = get(GravSlider, 'Value');
        set(GravEdit, 'String', num2str(num));
        
    end 
    function SnapCallback(varargin)
        
        F=getframe;
        [im,map]=frame2im(F);
        x=imshow(im);
        imsave(x);
    end

    function BrowseButtonPressed(varargin)
        filename = uigetfile({'*.dem';'*.*';'*.m'},'Select an input file');
        set(FileEdit, 'String', filename);
    end

    function StartButtonPressed(varargin)
                Val = get(ListBox, 'Value');
        
        switch Val
            case 1 % Rubber
                young = 0.1e7; % GPa
                density = 1.1e3; %kg/m3
            case 2 % Aluminium
                young = 69e7; % GPa
                density = 2.7e3; %kg/m3
            case 3 % Glass
                young = 50e7; % GPa
                density = 2.5e3; %kg/m3
            case 4 % Silicon
                young = 150e7; % GPa
                density = 2.33e3; %kg/m3
            case 5 % Steel
                young = 200e7; % GPa
                density = 7.82e3; %kg/m3
            case 6 % Polypropylene
                young = 1.5e7; % GPa
                density = 0.946e3; %kg/m3
        end

        ts = str2num(get(StepEdit,'String'));
        time = str2num(get(TimeEdit, 'String'));
        fric = str2num(get(FricEdit, 'String'));
        filename = get(FileEdit, 'String');

        gravity=str2num(get(GravEdit, 'String'));
        dem(filename,ts,time,gravity,young,density,fric);
    end

%% Simulator
function dem (file, timestep, time, gravity, young, density, mu)

global Partners Pindex np mp x_0 y_0 lx ly A gamma Pos Vel Acc Data...
    Force gk gm nx ny MassInertia pbcase;

A = 0.01;
gamma = 10;

plotupdate = 200;
plotcount=0;

Initialize(file);
drawnow
switch pbcase
    case 'on'
        for t=0:timestep:time
            Lattice(); 
            UpdatePart1(timestep); % Position Update and Velocity half step
            zero_Force(); % Setting the force vector to zero;
            make_Force(); % Calculating the interaction forces acting between two particles
            UpdatePart2(timestep); % Calculating the vel and acc of the particles
            periodic_bc(); % Check the boundary condition
            if (plotcount == plotupdate)
                Output();
                 set(CurrentTimeText, 'String', num2str(t));
                drawnow
                plotcount=0;
            end
            plotcount=plotcount+1;
        end
    case 'off'        
        for t=0:timestep:time
            Lattice(); 
            UpdatePart1(timestep); % Position Update and Velocity half step
            zero_Force(); % Setting the force vector to zero;
            make_Force(); % Calculating the interaction forces acting between two particles
            UpdatePart2(timestep); % Calculating the vel and acc of the particles
            if (plotcount == plotupdate)
                Output();
                 set(CurrentTimeText, 'String', num2str(t));
                drawnow
                plotcount=0;
            end
            plotcount=plotcount+1;
        end
end



    function Initialize (file)
        % file example
        %x_0,y_0,lx,ly
        % 3
        % posx,posy,velx,vely,omega,radius,spezies
        
        fid= fopen(file,'r');
        line=fgetl(fid);
        [values(1:4,1)]=sscanf(line,'%f, %f, %f, %f');
        x_0=values(1,1);
        y_0=values(2,1);
        lx=values(3,1);
        ly=values(4,1);
        line = fgetl(fid);
        np = sscanf(line,'%u');
        mp=0;
        Pos=zeros(np,2);
        Vel=zeros(np,3);
        for i=1:np
            line=fgetl(fid);
            [values(1:7,1)]=sscanf(line,'%f, %f, %f, %f, %f, %f, %f');
            Pos(i,1:2)=values(1:2,1);
            Vel(i,1:3)=values(3:5,1);
            Data(i,1)=values(6,1);
            Data(i,2)=values(7,1);
            Data(i,3)=4/3*pi*Data(i,1)^3*density;
            Data(i,4) = 2/5*Data(i,1)*Data(i,1)*Data(i,1);
            if Data(i,2)>0
                mp=mp+1;
            end
        end
        Acc=zeros(mp,3);
        MassInertia=[Data(1:mp,3) Data(1:mp,3) Data(1:mp,4)];
        line=fgetl(fid);
        pbcase=sscanf(line, '%s');
        fclose(fid);
        
        Rmin=Data(1,1);
        Rmax=Data(1,1);
        for i=1:np
            if(Data(i,1)<Rmin)
                Rmin=Data(i,1);
            end
            if (Data(i,1)>Rmax)
                Rmax=Data(i,1);
            end
        end
        gk=sqrt(2)*Rmin;
        gm=round(2*Rmax/gk)+1;
        nx=round(lx/gk)+1;
        ny=round(ly/gk)+1;
        Partners=zeros(np,30);
        Pindex=-ones(nx,ny);
        
        for i=1:np
            if  Data(i,2)== 0
                filledCircle(Pos(i,1:2),Data(i,1),50,'r');
            elseif Data(i,2)==1
                filledCircle(Pos(i,1:2),Data(i,1),50,'b');
            elseif Data(i,2)==2
                filledCircle(Pos(i,1:2),Data(i,1),50,'g');
            end
            hold on
        end
        set (gca,'XLIM',[x_0 x_0+lx],'YLIM',[y_0 y_0+ly],'Position',[0.35,0.075,0.6,0.8]);
        hold off
    end
    function Lattice()
        for i=1:np
            x=Pos(i,1);
            y=Pos(i,2);
            if(x>=x_0) && (x<x_0+lx) && (y>=y_0) && (y<y_0+ly)
                ix=round((x-x_0)/gk)+1;
                iy=round((y-y_0)/gk)+1;
                Pindex(ix,iy)=i;
                Partners(i,:)=0;
            end
        end
        for i=1:np
            x=Pos(i,1);
            y=Pos(i,2);
            if(x>=x_0) && (x<x_0+lx) && (y>=y_0) && (y<y_0+ly)
                ix=round((x-x_0)/gk)+1;
                iy=round((y-y_0)/gk)+1;
                for dx=-gm:gm
                    for dy=-gm:gm
                        k=Pindex(mod(ix+dx+nx,nx)+1,mod(iy+dy+ny,ny)+1);
                        if k>i
                            j=1;
                            while Partners(i,j)> 0
                                j=j+1;
                            end
                            Partners(i,j)=k;
                        end
                    end
                end
            end
        end
        Pindex=-ones(nx,ny);
    end
    function zero_Force()
       Force=zeros(mp,3);
    end
    function make_Force()
        for i=1:mp
            k=1;
            while k<30 && Partners(i,k)~=0
                pk=Partners(i,k);
                force(i,pk);
                k=k+1;
            end
        end
    end
    function force(P1,P2)
        dx=normalize(Pos(P1,1)-Pos(P2,1),lx);
        dy=normalize(Pos(P1,2)-Pos(P2,2),ly);
        rr=sqrt(dx*dx+dy*dy);
        r1=Data(P1,1);
        r2=Data(P2,1);
        xi=r1+r2-rr;
        if xi>0
            reff=(r1*r2)/(r1+r2);
            dvx=Vel(P1,1)-Vel(P2,1);
            dvy=Vel(P1,2)-Vel(P2,2);
            ex=dx/rr;
            ey=dy/rr;
            xidot=-(ex*dvx+ey*dvy);
            vtrel=-dvx*ey+dvy*ex+Vel(P1,3)*r1-Vel(P2,3)*r2;
            fn=sqrt(xi)*young*sqrt(reff)*(xi+A*xidot);
            ft=-gamma*vtrel;
            if(fn<0)
                fn=0;
            end
            if(ft<-mu*fn)
                ft=-mu*fn;
            end
            if(ft>mu*fn)
                ft=mu*fn;
            end
            if Data(P1,2)~=0
                Force(P1,:)=Force(P1,:)+[fn*ex-ft*ey,fn*ey+ft*ex,r1*ft];
            end
            if Data(P2,2)~=0
                Force(P2,:)=Force(P2,:)-[fn*ex-ft*ey,fn*ey+ft*ex,r2*ft];
            end
        end
    end
    function dx=normalize(dx,L)
        while(dx<-L/2)
            dx=dx+L;
        end
        while(dx>=L/2)
            dx=dx-L;
        end
    end
    function UpdatePart1(dt)
        %Updating Acceleration, Velocity and Position
        Pos(1:mp,:)=Pos(1:mp,:)+Vel(1:mp,1:2).*dt+1/2*Acc(1:mp,1:2).*dt.*dt;
        Vel(1:mp,:)=Vel(1:mp,:)+1/2*Acc.*dt;
    end
    function UpdatePart2(dt)
        Acc=Force./MassInertia;
        Acc(1:mp,2)=Acc(1:mp,2)-gravity;
        Vel(1:mp,:)=Vel(1:mp,:)+1/2*Acc.*dt;
    end
    function periodic_bc()
        for i=1:np
            while Pos(i,1)<x_0
                Pos(i,1)=Pos(i,1)+lx;
            end
            while Pos(i,1)>x_0+lx
                Pos(i,1)=Pos(i,1)-lx;
            end
            while Pos(i,2)<y_0
                Pos(i,2)=Pos(i,2)+ly;
            end
            while Pos(i,2)>y_0+ly
                Pos(i,2)=Pos(i,2)-ly;
            end
        end
    end
    function Output()
        %clf
        for i=1:np
            if Data(i,2)==0
                filledCircle(Pos(i,1:2),Data(i,1),50,'r');
            elseif Data(i,2)==1
                filledCircle(Pos(i,1:2),Data(i,1),50,'b');
            elseif Data(i,2)==2
                filledCircle(Pos(i,1:2),Data(i,1),50,'g');    
            end
            hold on
        end
        
        set (gca,'XLIM',[x_0 x_0+lx],'YLIM',[y_0 y_0+ly],'Position',[0.35,0.075,0.6,0.8]);
        hold off
    end
end
end