Nx=23;
Ny=132;
buffer.x= 0;
buffer.y= 10;
dx=0.00032609;
dy=0.00032609;

 [m n]=size(rER);
 cf.x = ceil((Nx - buffer.x*2-NPML(1)-NPML(2))/m); % Conversion factor to convert our real grid to our numerical grid
 cf.y = ceil((Ny - buffer.y*2-3-NPML(3)-NPML(4))/n);
 %Material Vectors
 ER = zeros([Nx Ny]);
 UR = zeros([Nx Ny]);
  
 % We Need to lay our real materials vectors over our numerical material
 % grid

%Fill our buffer regions
ER(1:buffer.x+NPML(1),:)=1;
ER(:,1:buffer.y+NPML(3)+2) = 1;
ER(Nx-buffer.x-NPML(2)+1:Nx,:) = 1;
ER(:,Ny-buffer.y-NPML(4)-1+1:Ny) = 1;

UR(1:buffer.x+NPML(1),:)=1;
UR(:,1:buffer.y+NPML(3)+2) = 1;
UR(Nx-buffer.x-NPML(2)+1:Nx,:) = 1;
UR(:,Ny-buffer.y-NPML(4)-1+1:Ny) = 1;


 
 % Lets place our real grid in proper location on numerical grid
 for x=0:m-1
    for y=0:n-1
      index.x = buffer.x+NPML(1) + x*cf.x+1;
      index.y = buffer.y+2+NPML(3) + y*cf.y+1;
      ER(index.x, index.y) = rER(x+1,y+1);
      UR(index.x, index.y) = rUR(x+1,y+1);
    end
 end
 
 %Fill our buffer regions
ER(1:buffer.x+NPML(1),:)=1;
ER(:,1:buffer.y+NPML(3)+2) = 1;
ER(Nx-buffer.x-NPML(2)+1:Nx,:) = 1;
ER(:,Ny-buffer.y-NPML(4)-1+1:Ny) = 1;

UR(1:buffer.x+NPML(1),:)=1;
UR(:,1:buffer.y+NPML(3)+2) = 1;
UR(Nx-buffer.x-NPML(2)+1:Nx,:) = 1;
UR(:,Ny-buffer.y-NPML(4)-1+1:Ny) = 1;

 
 for x=1 : cf.x : Nx
   for y=1 : Ny
    if(ER(x,y) == 0)
      ER(x,y) = ER(x,y-1);
    end
   
    if(UR(x,y) == 0)
      UR(x,y) = UR(x,y-1);
    end
   end
 end

 for y=1 : cf.y : Ny
  for x=1 : Nx
    if(ER(x,y) == 0)
      ER(x,y) = ER(x-1,y);
    end
   
    if(UR(x,y) == 0)
      UR(x,y) = UR(x-1,y);
    end
   end
 end
 
 
 
 %ER = ER';