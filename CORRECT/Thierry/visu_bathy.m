clear all ; close all
cd /Users/tpenduff/Documents/WORK/DATA/ORCA/ORCA_BATHY/BATHY_ORCA05/OPABAT05/CORRECT
c=load('/Users/tpenduff/pal_Testu.txt');
c=colormap(jet);

dir = '/Users/tpenduff/Documents/WORK/DATA/ORCA/ORCA_BATHY/BATHY_ORCA05/OPABAT05'; 

bat1=zeros(1442,1021);
bat2=zeros(1442,1021);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LECTURE 
namvarNC='Bathymetry';
%file = [ dir '/ORCA025_combined_etopo_gebco_coast_june2005.nc' ];
%file = [ dir '/ORCA025_combined_etopo_gebco_coast_corrected_oct05_G45_0.nc' ];
file = [ dir '/bathy_meter_treated_ORCA_R05_checked_coast_nobug.nc' ];
nc=netcdf(file,'nowrite');variables = var(nc);
i=1;lon=variables{i}(:);i=2;lat=variables{i}(:);[ny,nx]=size(lon);
bat1=nc{namvarNC}(:,:);


%file = [ dir '/ORCA025_combined_etopo_gebco_coast_corrected_oct05_G45_1.nc' ];
 file = [ dir '/bathy_meter_treated_ORCA_R05_checked_coast_nobug_tpenduff.nc' ];
nc=netcdf(file,'nowrite');variables = var(nc);
i=1;lon=variables{i}(:);i=2;lat=variables{i}(:);[ny,nx]=size(lon);
bat2=nc{namvarNC}(:,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LECTURE O25 FINALE PORU COMPARAISON
dir3= '/Users/tpenduff/Documents/WORK/DATA/ORCA/ORCA_BATHY/BATHY_ORCA025/'; 
file = [ dir3 '/ORCA025_combined_etopo_gebco_coast_corrected_oct05_G47.nc' ];
nc=netcdf(file,'nowrite');variables = var(nc);
i=1;lon3=variables{i}(:);i=2;lat3=variables{i}(:);[ny3,nx3]=size(lon3);
bat3=nc{namvarNC}(:,:);


%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% ALL PLOTS 3D
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

vx=[1:nx];vy=[1:ny]; 
figure(1);clf

nzones = 7;

for iplot=1:nzones

 if iplot==1;vx=[325:330];vy=[559 : 567]; V=[-1500 -150]  ; VIEW=[-63 64];  namzone='GIB' ;  end 
 if iplot==2;vx=[410:419];vy=[515 :525 ]; V=[-1800 -500]  ; VIEW=[-30 80];  namzone='DEN' ;  end 
 if iplot==3;vx=[396:403];vy=[547 :558 ]; V=[-1200  -400] ; VIEW=[-112 84]; namzone='FBC' ;  end 
 if iplot==4;vx=[265:275];vy=[485:498]  ; V=[-4700 -3800] ; VIEW=[-61 86];  namzone='VFZ' ;  end 
 if iplot==5;vy=[659:667];vx=[271:280]  ; V=[-1200 -000]  ; VIEW=[45 78];   namzone='BEM' ;  end 
 if iplot==6;vy=[524:535];vx=[8:25]     ; V=[-800  -000]  ; VIEW=[45 78];   namzone='WED' ;  end 
 if iplot==7;vx=[225:235]; vy=[135:145] ; V=[-20 0]       ; VIEW=[-63 64];  namzone='TOR' ;  end 

 if iplot==1;vx3=[650:657];vy3=[1118:1130];  end 
 if iplot==2;vx3=[820:838];vy3=[1030:1050];  end 
 if iplot==3;vx3=[792:805];vy3=[1095:1115];  end 
 if iplot==4;vx3=[530:550];vy3=[970:995]  ;  end 
 if iplot==5;vy3=[1318:1334];vx3=[542:558];  end 
 if iplot==6;vy3=[1046:1072];vx3=[15:50  ];  end 
 if iplot==7;vx3=[450:470];vy3=[270:290]  ;  end 

X=vy;Y=vx;A=-bat1(vx,vy);
[I J]=size(A);B=zeros(2*I,2*J);XX=zeros(1,2*J);YY=zeros(1,2*I);for i=1:I;for j=1:J;B (1+2*(i-1),1+2*(j-1))=A(i,j);B (2+2*(i-1),1+2*(j-1))=A(i,j);B (1+2*(i-1),2+2*(j-1))=A(i,j);B (2+2*(i-1),2+2*(j-1))=A(i,j);end;end
for j=1:J-1;XX(1+2*(j-1))=X(j);XX(2+2*(j-1))=X(j+1);end;XX(2*J-1)=X(J);XX(2*J)=X(J)+1;XX=XX-.5;
for i=1:I-1;YY(1+2*(i-1))=Y(i);YY(2+2*(i-1))=Y(i+1);end;YY(2*I-1)=Y(I);YY(2*I)=Y(I)+1;YY=YY-.5;
figure(1);subplot(nzones,3,1+3*(iplot-1)); surf(XX,YY,B);view(2);axis tight;shading flat;caxis(V);colormap(c);colorbar
ylabel('depth bat1')
figure(2+iplot-1);clf;subplot(221);   surf(XX,YY,B);caxis(V);colormap(c);colorbar;view(VIEW);axis tight;
ylabel('depth bat1')
AA=axis;


X=vy;Y=vx;A=-bat2(vx,vy);
[I J]=size(A);B=zeros(2*I,2*J);XX=zeros(1,2*J);YY=zeros(1,2*I);for i=1:I;for j=1:J;B (1+2*(i-1),1+2*(j-1))=A(i,j);B (2+2*(i-1),1+2*(j-1))=A(i,j);B (1+2*(i-1),2+2*(j-1))=A(i,j);B (2+2*(i-1),2+2*(j-1))=A(i,j);end;end
for j=1:J-1;XX(1+2*(j-1))=X(j);XX(2+2*(j-1))=X(j+1);end;XX(2*J-1)=X(J);XX(2*J)=X(J)+1;XX=XX-.5;
for i=1:I-1;YY(1+2*(i-1))=Y(i);YY(2+2*(i-1))=Y(i+1);end;YY(2*I-1)=Y(I);YY(2*I)=Y(I)+1;YY=YY-.5;
figure(1);subplot(nzones,3,2+3*(iplot-1)); surf(XX,YY,B);view(2);axis tight;shading flat;caxis(V);colormap(c);colorbar
figure(1);subplot(nzones,3,3+3*(iplot-1));contour(vy,vx,-bat2(vx,vy),[V(1):200:V(2)]);hold on;grid
clear dif;dif=bat2(vx,vy)-bat1(vx,vy);dif(dif==0)=NaN;pcolor(vy-.5,vx-.5,dif);shading flat;colormap(c);colorbar;
;axis tight;axis(AA);
shading flat;ylabel('depth bat2-bat1')

figure(2+iplot-1);subplot(222);      surf(XX,YY,B);caxis(V);colormap(c);colorbar;view(VIEW);axis tight
ylabel('depth bat2')

X3=vy3;Y3=vx3;A3=-bat3(vx3,vy3);
[I3 J3]=size(A3);B3=zeros(2*I3,2*J3);XX3=zeros(1,2*J3);YY3=zeros(1,2*I3);for i=1:I3;for j=1:J3;B3(1+2*(i-1),1+2*(j-1))=A3(i,j);B3(2+2*(i-1),1+2*(j-1))=A3(i,j);B3(1+2*(i-1),2+2*(j-1))=A3(i,j);B3(2+2*(i-1),2+2*(j-1))=A3(i,j);end;end
for j=1:J3-1;XX3(1+2*(j-1))=X3(j);XX3(2+2*(j-1))=X3(j+1);end;XX3(2*J3-1)=X3(J3);XX3(2*J3)=X3(J3)+1;XX3=XX3-.5;
for i=1:I3-1;YY3(1+2*(i-1))=Y3(i);YY3(2+2*(i-1))=Y3(i+1);end;YY3(2*I3-1)=Y3(I3);YY3(2*I3)=Y3(I3)+1;YY3=YY3-.5;
figure(2+iplot-1);subplot(224);      surf(XX3,YY3,B3);caxis(V);colormap(c);colorbar;view(VIEW);axis tight
AA3=axis;axis([AA3(1) AA3(2) AA3(3) AA3(4) AA(5) AA(6)]);
ylabel('depth bat3 (025 G47)')

%figure(6);clf;surf(lon3(YY3-.5,XX3-.5),lat3(YY3-.5,XX3-.5),B3);caxis(V);colormap(c);colorbar;view(VIEW);axis tight
end




%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% ALL PLOTS 3D COLORES PAR k
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
d=load('../depw.txt');
vx=[1:nx];vy=[1:ny]; 

for iplot=1:nzones

 if iplot==1;vx=[325:330];vy=[559 : 567]; V=[-1500 -150]  ; VIEW=[-63 64];  namzone='GIB' ;  end 
 if iplot==2;vx=[410:419];vy=[515 :525 ]; V=[-1800 -500]  ; VIEW=[-30 80];  namzone='DEN' ;  end 
 if iplot==3;vx=[396:403];vy=[547 :558 ]; V=[-1200  -400] ; VIEW=[-112 84]; namzone='FBC' ;  end 
 if iplot==4;vx=[265:275];vy=[485:498]  ; V=[-4700 -3800] ; VIEW=[-61 86];  namzone='VFZ' ;  end 
 if iplot==5;vy=[659:667];vx=[271:280]  ; V=[-1200 -000]  ; VIEW=[45 78];   namzone='BEM' ;  end 
 if iplot==6;vy=[524:535];vx=[8:25]     ; V=[-800  -000]  ; VIEW=[45 78];   namzone='WED' ;  end 
 if iplot==7;vx=[225:235]; vy=[135:145] ; V=[-20 0]       ; VIEW=[-63 64];  namzone='TOR' ;  end 

 if iplot==1;vx3=[650:657];vy3=[1118:1130];  end 
 if iplot==2;vx3=[820:838];vy3=[1030:1050];  end 
 if iplot==3;vx3=[792:805];vy3=[1095:1115];  end 
 if iplot==4;vx3=[530:550];vy3=[970:995]  ;  end 
 if iplot==5;vy3=[1318:1334];vx3=[542:558];  end 
 if iplot==6;vy3=[1046:1072];vx3=[15:50  ];  end 
 if iplot==7;vx3=[450:470];vy3=[270:290]  ;  end 

figure(nzones+1+iplot);clf
X=vy;Y=vx;A=-bat2(vx,vy);
[I J]=size(A);B=zeros(2*I,2*J);XX=zeros(1,2*J);YY=zeros(1,2*I);for i=1:I;for j=1:J;B (1+2*(i-1),1+2*(j-1))=A(i,j);B (2+2*(i-1),1+2*(j-1))=A(i,j);B (1+2*(i-1),2+2*(j-1))=A(i,j);B (2+2*(i-1),2+2*(j-1))=A(i,j);end;end
for j=1:J-1;XX(1+2*(j-1))=X(j);XX(2+2*(j-1))=X(j+1);end;XX(2*J-1)=X(J);XX(2*J)=X(J)+1;XX=XX-.5;
for i=1:I-1;YY(1+2*(i-1))=Y(i);YY(2+2*(i-1))=Y(i+1);end;YY(2*I-1)=Y(I);YY(2*I)=Y(I)+1;YY=YY-.5;
subplot(221);       surf(XX,YY,B);caxis(V);colormap(c);colorbar;view(VIEW);axis tight; title('depth bat2')

%%% Calcul de kabove = indice de niveau k juste au-dessus topo Pstep locale.
%%% Calcul de AFS    = profondeur Fstep   juste au dessus de la topo Pstep locale.
X=vy;Y=vx;A=-bat2(vx,vy);
[I J]=size(A);kabove=zeros(I,J);AFS=zeros(I,J);
%for i=1:I;for j=1:J;if -A(i,j)>0;kabove(i,j)=max(find(-A(i,j)./d(2:46)>=1))+2;end;if kabove(i,j)>0;AFS(i,j)=-d(kabove(i,j));end;end;end;
 for i=1:I;for j=1:J;
if -A(i,j)>0;
 %kabove(i,j)= min(find(-A(i,j)./d(2:46)<1))+1;
  kabove(i,j)= min(find(-A(i,j)./d(2:46)<1));
 if kabove(i,j)>0;AFS(i,j)=-d(kabove(i,j));end;
end;
 end;end;

X=vy;Y=vx;A=-bat2(vx,vy);
[I J]=size(A);B=zeros(2*I,2*J);XX=zeros(1,2*J);YY=zeros(1,2*I);for i=1:I;for j=1:J;B (1+2*(i-1),1+2*(j-1))=A(i,j);B (2+2*(i-1),1+2*(j-1))=A(i,j);B (1+2*(i-1),2+2*(j-1))=A(i,j);B (2+2*(i-1),2+2*(j-1))=A(i,j);end;end
B0=B;
X=vy;Y=vx;A=kabove;
[I J]=size(A);B=zeros(2*I,2*J);XX=zeros(1,2*J);YY=zeros(1,2*I);for i=1:I;for j=1:J;B (1+2*(i-1),1+2*(j-1))=A(i,j);B (2+2*(i-1),1+2*(j-1))=A(i,j);B (1+2*(i-1),2+2*(j-1))=A(i,j);B (2+2*(i-1),2+2*(j-1))=A(i,j);end;end
B1=B;
X=vy;Y=vx;A=AFS;
[I J]=size(A);B=zeros(2*I,2*J);XX=zeros(1,2*J);YY=zeros(1,2*I);for i=1:I;for j=1:J;B (1+2*(i-1),1+2*(j-1))=A(i,j);B (2+2*(i-1),1+2*(j-1))=A(i,j);B (1+2*(i-1),2+2*(j-1))=A(i,j);B (2+2*(i-1),2+2*(j-1))=A(i,j);end;end
B2=B;
for j=1:J-1;XX(1+2*(j-1))=X(j);XX(2+2*(j-1))=X(j+1);end;XX(2*J-1)=X(J);XX(2*J)=X(J)+1;XX=XX-.5;
for i=1:I-1;YY(1+2*(i-1))=Y(i);YY(2+2*(i-1))=Y(i+1);end;YY(2*I-1)=Y(I);YY(2*I)=Y(I)+1;YY=YY-.5;

subplot(222);surf(XX,YY,B0,-B1);view(VIEW);axis tight;shading faceted;colormap(c);colorbar  ;title('k index below h_{Pstep} along bat2')
subplot(223);surf(XX,YY,B2,B2) ;view(VIEW);axis tight;shading faceted;colormap(c);caxis(V);colorbar;title('h_{Fstep} below h_{Pstep}');view(2)
subplot(224);surf(XX,YY,B0,B0-B2) ;view(VIEW);axis tight;shading faceted;colormap(c);colorbar;title('h_{Pstep} - h_{Pstep}');view(2)

end

figure(15);clf;surf(-bat2);shading flat;view(2);colorbar;axis tight;CA=caxis;
figure(14);clf;surf(-bat1);shading flat;view(2);caxis(CA);colorbar;axis tight;

tmp=bat1-bat2;tmp(tmp>5000)=NaN;%tmp(tmp==0)=NaN;[t1 t2]=size(tmp);
figure(16);clf;surf(lon,lat,abs(tmp));shading flat;view(3);caxis([-1 1]);colorbar;axis tight;hold on;contour(lon,lat,bat2, [1 1],'k');
view(-30,30);axis([-180 180 -77 89.885 0 1000])

figure(17);clf;hist(reshape(tmp,1,t1*t2),50)

% figure(1) ;orient landscape;print -dpsc PLOTS/ALL.ps
% figure(2) ;orient landscape;print -dpsc PLOTS/GIB1.ps
% figure(3) ;orient landscape;print -dpsc PLOTS/DEN1.ps
% figure(4) ;orient landscape;print -dpsc PLOTS/FBC1.ps
% figure(5) ;orient landscape;print -dpsc PLOTS/VFZ1.ps
% figure(6) ;orient landscape;print -dpsc PLOTS/BEM1.ps
% figure(7) ;orient landscape;print -dpsc PLOTS/WED1.ps
% figure(8) ;orient landscape;print -dpsc PLOTS/GIB2.ps
% figure(9) ;orient landscape;print -dpsc PLOTS/DEN2.ps
% figure(10);orient landscape;print -dpsc PLOTS/FBC2.ps
% figure(11);orient landscape;print -dpsc PLOTS/VFZ2.ps
% figure(12);orient landscape;print -dpsc PLOTS/BEM2.ps
% figure(13);orient landscape;print -dpsc PLOTS/WED2.ps







%%%%%%%%%%%%%%%%%%%%% RUBBISH
%subplot(2,2,1);surf(A)           ;view(52,38)
%subplot(2,2,2);surf(-kabove)     ;view(52,38)
%subplot(2,2,3);surf(AFS)         ;view(52,38)
%subplot(2,2,4);surf(A-AFS)       ;view(52,38)








%figure(1);clf;hold on
%i=0;
%for zz=1:20:1000
%i=i+1;
%zzz=d(max(find(d/zz<=1))+1)
%plot(-zz,-zzz,'.')
%end
%plot(-zzz,-d(1:20),'r.')

%[I J]=size(A);B=zeros(2*I,2*J);XX=zeros(1,2*J);YY=zeros(1,2*I);for i=1:I;for j=1:J;B (1+2*(i-1),1+2*(j-1))=A(i,j);B (2+2*(i-1),1+2*(j-1))=A(i,j);B (1+2*(i-1),2+2*(j-1))=A(i,j);B (2+2*(i-1),2+2*(j-1))=A(i,j);end;end
%for j=1:J-1;XX(1+2*(j-1))=X(j);XX(2+2*(j-1))=X(j+1);end;XX(2*J-1)=X(J);XX(2*J)=X(J)+1;XX=XX-.5;
%for i=1:I-1;YY(1+2*(i-1))=Y(i);YY(2+2*(i-1))=Y(i+1);end;YY(2*I-1)=Y(I);YY(2*I)=Y(I)+1;YY=YY-.5;
%X=vy;Y=vx;A=-bat1(vx,vy);
%[I J]=size(A);B=zeros(2*I,2*J);XX=zeros(1,2*J);YY=zeros(1,2*I);for i=1:I;for j=1:J;B (1+2*(i-1),1+2*(j-1))=A(i,j);B (2+2*(i-1),1+2*(j-1))=A(i,j);B (1+2*(i-1),2+2*(j-1))=A(i,j);B (2+2*(i-1),2+2*(j-1))=A(i,j);end;end
%for j=1:J-1;XX(1+2*(j-1))=X(j);XX(2+2*(j-1))=X(j+1);end;XX(2*J-1)=X(J);XX(2*J)=X(J)+1;XX=XX-.5;
%for i=1:I-1;YY(1+2*(i-1))=Y(i);YY(2+2*(i-1))=Y(i+1);end;YY(2*I-1)=Y(I);YY(2*I)=Y(I)+1;YY=YY-.5;


