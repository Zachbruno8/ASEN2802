clear;
clc;
close all;

%% inputs

N = 4;
Omega0 = 3; % lb/in
Omega1 = 3; % lb/in
LengthUnit = "in";
ForceUnit = "lb";
L = 10; % in
MofElasticity = 10E7;% lb/in^2
Width = .5; % in
Thicness = .12; % in

%% Calculating

fx = @(x) Omega0 + x * (Omega1-Omega0) / L; % lift equation
gx = @(x) x ^ 2 * Omega0 / 2 + x ^ 3 * (Omega1-Omega0) / 3 / L; % integral of lift equation * x

Deltax = L / N;

figure('units','inch','position',[1,1,13,8])
hold on
axis off
% Plotting lines
plot([0 L],[Omega0 Omega1],'r') % plotting lift line
plot([0 0],[0 Omega0],'b--') % Plotting leftmost vertical line
plot([0 L],[0 0],'k') % bottom line at y = 0
% window limits
ylim([-.5*Omega0 1.2*Omega0]) % Setting y limits
xlim([-.2*L 1.2*L]) % setting x limits 
% inserting text
text(L*-.02,0,append('L= ',string(L),' ',LengthUnit),'HorizontalAlignment','right') % Length
text(0,Omega0*1.1,append('\omega_0 = ',string(Omega0),' ',ForceUnit,'/',LengthUnit)) % Omega 0
text(L*1.01,Omega1,append('\omega_1 = ',string(Omega1),' ',ForceUnit,'/',LengthUnit)) % Omega 1


%% For loop to calculate each area and position

% Preallocating
FR = zeros(1,N);
x = zeros(1,N+1);
MassCenter = zeros(1,N);

for i = 1:N
    x(i+1) = x(i) + Deltax; % Adding delta x for next loop
    a = fx(x(i)); % left height
    b = fx(x(i)+Deltax);% right height
    FR(i) = .5 * (a + b) * Deltax; % area of trapezoid equation
    MassCenter(i) = (gx(x(i+1))-gx(x(i)))/FR(i); % center of mass
    plot([MassCenter(i) MassCenter(i)],[0 -.2*Omega0],'r',Marker="v") % Plotting force arrows
    text(MassCenter(i),-.27*Omega0,append('F',string(i),'= ',string(FR(i)),' ',ForceUnit,newline,'x',string(i),'= ',string(MassCenter(i)),' ',string(LengthUnit)),"HorizontalAlignment","center") % Displaying force and location
    plot([x(i+1) x(i+1)],[0 fx(x(i+1))],'b--') % Plotting vertical lines
end

% Finding I
I = Width * Thicness^3 /12;
% Finding Feq and xeq
Feq = sum(FR);
xeq = sum(MassCenter .* FR) / Feq;
v = (Feq * xeq^2)/(6 * MofElasticity * I) * (3 * L - xeq);

% Printing results
fprintf("Feq = %g %s at xeq = %g %s\n",Feq,ForceUnit,xeq,LengthUnit)
fprintf("   FR%g= %g %s at x%g= %g %s (%g)\n",[1:N;FR;repmat(ForceUnit,1,N);1:N;MassCenter;repmat(LengthUnit,1,N);MassCenter-x(1:N)])
fprintf("\n v max is %g %s\n",v,LengthUnit)

if N==4
%% Finding wiffle tree values
RodLength = zeros(1,3);
RodLength(2) = MassCenter(2) - MassCenter(1);
RodLength(3) = MassCenter(4) - MassCenter(3);
c = RodLength(2)*FR(2)/(FR(2)+FR(1));
d = RodLength(2)-c;
e = RodLength(3)*FR(4)/(FR(4)+FR(3));
f = RodLength(3)-e;
a = xeq - c - MassCenter(1);
b = MassCenter(4)- f - xeq;
RodLength(1) = a + b;
% Displaying
fprintf("\na = %g\nb = %g\nc = %g\nd = %g\ne = %g\nf = %g\n\n",a,b,c,d,e,f)
fprintf("Rod %g = %g\n",[1 2 3;RodLength])

%% Plotting
figure('units','inch','position',[1,1,13,8])
axis off
hold on
ylim([-L/10 L/10]) % Setting y limits
xlim([-.2*L 1.2*L]) % setting x limits 
% wing
plot([0 L],[L/20 L/20],"k","LineWidth",5)
% Bars 1 2 and 3
plot([MassCenter(1) MassCenter(2)],[0 0],"b","LineWidth",20)
plot([MassCenter(3) MassCenter(4)],[0 0],"b","LineWidth",20)
plot([MassCenter(1)+c MassCenter(3)+e],[-L/20 -L/20],"b","LineWidth",20)
% F_12 and F_13 lines
plot([MassCenter(1)+c MassCenter(1)+c],[L/100 -L/20],"k")
plot([MassCenter(3)+e MassCenter(3)+e],[L/100 -L/20],"k")
% F_R lines
plot([MassCenter(1) MassCenter(1)],[0 L/20],"k")
plot([MassCenter(2) MassCenter(2)],[0 L/20],"k")
plot([MassCenter(3) MassCenter(3)],[0 L/20],"k")
plot([MassCenter(4) MassCenter(4)],[0 L/20],"k")
plot([xeq xeq],[-L/25 -L/15],"k")
% adding text
% bar 2
text(MassCenter(1)+.5*c,L/100,append("c= ",newline,string(c)),"HorizontalAlignment","center","VerticalAlignment","baseline")
text(MassCenter(2)-.5*d,L/100,append("d= ",newline,string(d)),"HorizontalAlignment","center","VerticalAlignment","baseline")
% bar 3
text(MassCenter(3)+.5*e,L/100,append("e= ",newline,string(e)),"HorizontalAlignment","center","VerticalAlignment","baseline")
text(MassCenter(4)-.5*f,L/100,append("f= ",newline,string(f)),"HorizontalAlignment","center","VerticalAlignment","baseline")
% bar 1
text(MassCenter(1)+c+.5*a,-9*L/200,append("a= ",string(a)),"HorizontalAlignment","center","VerticalAlignment","baseline")
text(MassCenter(3)+e-.5*b,-9*L/200,append("b= ",string(b)),"HorizontalAlignment","center","VerticalAlignment","baseline")
% bars
text(MassCenter(2)+L/200,0,append("Rod 2= ",string(RodLength(2))),"HorizontalAlignment","left","VerticalAlignment","middle")
text(MassCenter(4)+L/200,0,append("Rod 3= ",string(RodLength(3))),"HorizontalAlignment","left","VerticalAlignment","middle")
text(MassCenter(3)+e+L/200,-L/20,append("Rod 1= ",string(RodLength(1))),"HorizontalAlignment","left","VerticalAlignment","middle")
% F_R
text(xeq,-L/15,append("Fr= ",string(Feq),newline,"at x= ",string(xeq)),"HorizontalAlignment","center","VerticalAlignment","top")
text(MassCenter(1),L/22,append(" x1= ",string(MassCenter(1))),"HorizontalAlignment","left","VerticalAlignment","middle")
text(MassCenter(2),L/22,append(" x2= ",string(MassCenter(2))),"HorizontalAlignment","left","VerticalAlignment","middle")
text(MassCenter(3),L/22,append(" x3= ",string(MassCenter(3))),"HorizontalAlignment","left","VerticalAlignment","middle")
text(MassCenter(4),L/22,append(" x4= ",string(MassCenter(4))),"HorizontalAlignment","left","VerticalAlignment","middle")
print(gcf,'Wiffle_Tree.png','-dpng','-r150')
end

% vx = @(x) (Omega0 + ((Omega1-Omega0)/L) * L ^2 /2)*((Omega0+2*Omega1)/(3*(Omega0+Omega1))*L)*x^2/2;
vx = @(x) (Omega0 + ((Omega1-Omega0)/L) * L ^2 /2)*((Omega0+2*Omega1)/(3*(Omega0+Omega1))*L)*x.^2/2;
if Omega0 == Omega1
Mx = @(x) Omega0 * x .^ 2 / 2; % lift equation
Sigmax = @(x) -Omega0 * x .^ 2 * Thicness / 4 / I;
Epsilonx = @(x) -Omega0 * x .^ 2 * Thicness / 4 / I / MofElasticity;
Vx = @(x) Omega0 * x .^ 4 / 24 / MofElasticity / I;
Vx1 = @(x) Feq/6/MofElasticity/I*(-x.^3+3*xeq*x.^2);
Vx2 = @(x) Feq/(6*MofElasticity*I)*(3*xeq^2*x-xeq^3);
v = @(x) Omega0 * x .^ 5 / 120 / MofElasticity / I;
Lx = linspace(0,L,100);

figure
hold on
title("Deflection")
% plot(Lx,Vx(Lx))
% plot(Lx,v(Lx))
plot(Lx(1,1:50),Vx1(Lx(1,1:50)),"b")
plot(Lx(1,50:100),Vx2(Lx(1,50:100)),"b")
plot(Lx,vx(Lx))
ylabel("v (in)")
xlabel("L (in)")
print(gcf,'Deflection.png','-dpng','-r150')

figure
hold on
title("Stress")
plot(Lx,Sigmax(Lx)-Sigmax(L))
ylabel("sigma (lb/in^2)")
xlabel("L (in)")
print(gcf,'Stress.png','-dpng','-r150')

figure
hold on
title("Strain")
plot(Lx,Epsilonx(Lx)-Epsilonx(L))
ylabel("Strain")
xlabel("L (in)")
print(gcf,'Strain.png','-dpng','-r150')

end