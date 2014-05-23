%% PCA Validation
% This is a reference script for how to obtain basis vectors and scores under different conditions when using SVD. 
% Use it to guide you through the forest when darkness descends. With abundant graphical representations put your mind at ease.
%% Initialization
clear;clc;
%% Create Ellipsoid
[x y z]=ellipsoid(0,0,0,15,10,5,20);
%% Rotation
ax=random('unif',0,pi);
ay=random('unif',0,pi);
az=random('unif',0,pi);
% First Rotation - About X Axis
x1=x;
y1=y*cos(ax)-z*sin(ax);
z1=y*sin(ax)+z*cos(ax);
% First Rotation - About X Axis
x2=x1*cos(ay)+z1*sin(ay);
y2=y1;
z2=-x1*sin(ay)+z1*cos(ay);
% First Rotation - About X Axis
x3=x2*cos(az)-y2*sin(az);
y3=x2*sin(az)+y2*cos(az);
z3=z2;
%% Translation
t1=random('unif',10,20);
t2=random('unif',10,20);
t3=random('unif',10,20);
x4=x3+t1;
y4=y3+t2;
z4=z3+t3;
%% Plot
figure(1);
line([0 0],[0 0],[0 30],'color','black','LineWidth',3);hold on;
line([0 0],[0 30],[0 0],'color','black','LineWidth',3);
line([0 30],[0 0],[0 0],'color','black','LineWidth',3);
text(0,0,32,'Z');text(32,0,0,'X');text(0,32,0,'Y')
surf(x4,y4,z4);axis equal;colormap copper;axis off;
%% Vectorize
X(1:21*21)=x4(1:21*21);
Y(1:21*21)=y4(1:21*21);
Z(1:21*21)=z4(1:21*21);
list=[X',Y',Z'];
%% Column Center - Parameter Center
M=repmat(mean(list,1),441,1);
Data=list-M;
%% PCA - Princomp
[C K L]=princomp(list);
%Plot
figure(2);
line([0 0],[0 0],[0 30],'color','black','LineWidth',3);hold on;
line([0 0],[0 30],[0 0],'color','black','LineWidth',3);
line([0 30],[0 0],[0 0],'color','black','LineWidth',3);
text(0,0,32,'Z');text(32,0,0,'X');text(0,32,0,'Y')
surf(x4,y4,z4);axis equal;colormap copper;axis off;
line([t1 t1+20*C(1,1)],[t2 t2+20*C(2,1)],[t3 t3+20*C(3,1)],'color','red','LineWidth',3);
line([t1 t1+15*C(1,2)],[t2 t2+15*C(2,2)],[t3 t3+15*C(3,2)],'color','red','LineWidth',3);
line([t1 t1+10*C(1,3)],[t2 t2+10*C(2,3)],[t3 t3+10*C(3,3)],'color','red','LineWidth',3);
text(t1+20*C(1,1),t2+20*C(2,1),t3+20*C(3,1),'C1');
text(t1+15*C(1,2),t2+15*C(2,2),t3+15*C(3,2),'C2');
text(t1+10*C(1,3),t2+10*C(2,3),t3+10*C(3,3),'C3')
%% PCA - SVD
[U S V]=svd(Data,0);
%Plot
figure(3);
line([0 0],[0 0],[0 30],'color','black','LineWidth',3);hold on;
line([0 0],[0 30],[0 0],'color','black','LineWidth',3);
line([0 30],[0 0],[0 0],'color','black','LineWidth',3);
text(0,0,32,'Z');text(32,0,0,'X');text(0,32,0,'Y')
surf(x4,y4,z4);axis equal;colormap copper;axis off;
line([t1 t1+20*V(1,1)],[t2 t2+20*V(2,1)],[t3 t3+20*V(3,1)],'color','red','LineWidth',3);
line([t1 t1+15*V(1,2)],[t2 t2+15*V(2,2)],[t3 t3+15*V(3,2)],'color','red','LineWidth',3);
line([t1 t1+10*V(1,3)],[t2 t2+10*V(2,3)],[t3 t3+10*V(3,3)],'color','red','LineWidth',3);
text(t1+20*V(1,1),t2+20*V(2,1),t3+20*V(3,1),'C1');
text(t1+15*V(1,2),t2+15*V(2,2),t3+15*V(3,2),'C2');
text(t1+10*V(1,3),t2+10*V(2,3),t3+10*V(3,3),'C3');
%% PCA - PCA Lanczos Method (Code from Dave)
[U1 S1 V1]=pca(Data,3);
%Plot
figure(4);
line([0 0],[0 0],[0 30],'color','black','LineWidth',3);hold on;
line([0 0],[0 30],[0 0],'color','black','LineWidth',3);
line([0 30],[0 0],[0 0],'color','black','LineWidth',3);
text(0,0,32,'Z');text(32,0,0,'X');text(0,32,0,'Y')
surf(x4,y4,z4);axis equal;colormap copper;axis off;
line([t1 t1+20*V1(1,1)],[t2 t2+20*V1(2,1)],[t3 t3+20*V1(3,1)],'color','red','LineWidth',3);
line([t1 t1+15*V1(1,2)],[t2 t2+15*V1(2,2)],[t3 t3+15*V1(3,2)],'color','red','LineWidth',3);
line([t1 t1+10*V1(1,3)],[t2 t2+10*V1(2,3)],[t3 t3+10*V1(3,3)],'color','red','LineWidth',3);
text(t1+20*V1(1,1),t2+20*V1(2,1),t3+20*V1(3,1),'C1');
text(t1+15*V1(1,2),t2+15*V1(2,2),t3+15*V1(3,2),'C2');
text(t1+10*V1(1,3),t2+10*V1(2,3),t3+10*V1(3,3),'C3');
%% PCA - PCA Lanczos Method (Code from Dave as Dave Does It)
W=list';
Wm=W-repmat(mean(W, 2), [1 size(W,2)]);
[U2 S2 V2]=pca(Wm,3);
%Plot
figure(5);
line([0 0],[0 0],[0 30],'color','black','LineWidth',3);hold on;
line([0 0],[0 30],[0 0],'color','black','LineWidth',3);
line([0 30],[0 0],[0 0],'color','black','LineWidth',3);
text(0,0,32,'Z');text(32,0,0,'X');text(0,32,0,'Y')
surf(x4,y4,z4);axis equal;colormap copper;axis off;
line([t1 t1+20*U2(1,1)],[t2 t2+20*U2(2,1)],[t3 t3+20*U2(3,1)],'color','red','LineWidth',3);
line([t1 t1+15*U2(1,2)],[t2 t2+15*U2(2,2)],[t3 t3+15*U2(3,2)],'color','red','LineWidth',3);
line([t1 t1+10*U2(1,3)],[t2 t2+10*U2(2,3)],[t3 t3+10*U2(3,3)],'color','red','LineWidth',3);
text(t1+20*U2(1,1),t2+20*U2(2,1),t3+20*U2(3,1),'C1');
text(t1+15*U2(1,2),t2+15*U2(2,2),t3+15*U2(3,2),'C2');
text(t1+10*U2(1,3),t2+10*U2(2,3),t3+10*U2(3,3),'C3');