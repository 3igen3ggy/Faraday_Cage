clear
xdel(winsid())
sec_kind=%t
ani=%t
//Parameters
dx=50;//domain width (dim of a cage in x direction)
dy=50;//domain height (dim of a cage in y direction)
nx=201;
ny=201;
x=linspace(0,dx,nx);
y=linspace(0,dy,ny);
V=zeros(ndgrid(x,y))
Va=zeros(ndgrid(x,y))
V2=0  //right
V1=40    //left
V3=20  //up
V4=20 //down
Va(1,:)=V1
Va(nx,:)=V2
Va(:,1)=V3
Va(:,ny)=V4
flag=%t
acc=10^(-8)
i=1

k=1000;
Ax = 2*dx/k;
Ay = 2*dy/k;
function res = fbc2(X)
  res = X;
  res(:,1) = (4*X(:,2)-X(:,3)+Ay*V3)/(3+Ay);
  res(:,$) = (4*X(:,$-1)-X(:,$-2)+Ay*V4)/(3+Ay);
endfunction

number_of_rods=128; //4,8,16,32
space_between=0; //not to much, maximum depends on other parameters
temp=floor((nx-space_between*(number_of_rods-1)-number_of_rods)/2);
for j=1:(number_of_rods/4)
    rods(4*j-3,:)=[temp,temp+space_between*(j-1)+j-1];
    rods(4*j-2,:)=[temp+(space_between+1)*number_of_rods/4-space_between*(j-1)-j+1,temp];
    rods(4*j,:)=[temp+space_between*(j-1)+j-1,temp+(space_between+1)*number_of_rods/4];
    rods(4*j-1,:)=[temp+(space_between+1)*number_of_rods/4,temp+(space_between+1)*number_of_rods/4-space_between*(j-1)-j+1];
end

for (j=1:number_of_rods)
    R(j)=(rods(j,1)+ny*(rods(j,2)-1));
    end
    if ani then
    a=scf;
    a.color_map=jetcolormap(64)
    grayplot(x,y,V)
    isoview on
    f=gca();
    j=0;
end

while flag do
    Va(2:$-1,2:$-1)=0.25*(V(1:$-2,2:$-1)+V(3:$,2:$-1)+V(2:$-1,1:$-2)+V(2:$-1,3:$))
    Va(R)=0
    if sec_kind then 
        Va=fbc2(Va)
    end
    dV=max(abs((Va-V)./V))
    V=Va
    i=i+1
    if(ani) then
        j=j+1;
        if j==20 then
            f.children.data.z=V;
            j=0
        end
    end
    if dV<acc || i==10^10 then 
        flag=%f end   
    end

if ~ani then
    a=scf;
    a.color_map=jetcolormap(64)
    grayplot(x,y,V)
    isoview on
end

gradx = zeros(V);
grady = gradx;
gradx(2:$-1,2:$-1) = ( V(3:$,2:$-1) - V(1:$-2,2:$-1) ) /2/dx;
grady(2:$-1,2:$-1) = ( V(2:$-1,3:$) - V(2:$-1,1:$-2) ) /2/dy;
gradx(1,:) = ( V(2,:) - V(1,:) )/dx;
gradx($,:) = ( V($,:) - V($-1,:) )/dx;
grady(:,1) = ( V(:,2) - V(:,1) )/dy;
grady(:,$) = ( V(:,$) - V(:,$-1) )/dy;
grady(1,2:$-1) = ( V(1,3:$) - V(1,1:$-2) )/dx;
grady($,2:$-1) = ( V($,3:$) - V($,1:$-2) )/dx;
gradx(2:$-1,1) = ( V(3:$,1) - V(1:$-2,1) )/dy;
gradx(2:$-1,$) = ( V(3:$,$) - V(1:$-2,$) )/dy;

every = 5;

L = sqrt(gradx.^2+grady.^2);
Lmax = 0.01;
coef = L/Lmax;
coef(coef<1) = 1;
ngx = gradx./coef;
ngy = grady./coef;
scf();
champ(x(1:every:$),y(1:every:$),-ngx(1:every:$,1:every:$),-ngy(1:every:$,1:every:$))
ee = gce();
ee.arrow_size = 0.8;
xtitle('Clipped arrows','x','y')
isoview on;
