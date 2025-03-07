function [count,PWxList,PWyList,signList]=pw_finder_withsign(z)
% pinwheel finder optimized for theoretical orientation domains



% output variables
count=0;
PWxList=[];
PWyList=[];
signList=[];

% for avoidance of calc. pinwheel twice
doppeltestx=0; 
doppeltesty=0;
doppeldet=0;


% calculate contour lines
[sy,sx]  =size(z);
[grx,gry]=gradient(real(z));
[gix,giy]=gradient(imag(z));

aa=contourc(real(z),[0 0]);
bb=contourc(imag(z),[0 0]);
if ~isempty(aa) && ~isempty(bb)
    chre=getcontourlines(aa);
    lre=length(chre);
    chim = getcontourlines(bb);
    lim=length(chim);
else
    return
end

%% calculate distances & rectangles in advance!
Xim=cell(1,lim);
Yim=cell(1,lim);
Nim=zeros(1,lim);
Dist2=zeros(1,lim);
im_4=zeros(lim,4);
 for s= 1:lim

     xim =chim(s).x;
     yim =chim(s).y;
     l=length(xim);

     Xim{s} =xim(1:l);
     Yim{s} =yim(1:l);  %'NaN abschneiden)
     Nim(s)     = l;

     if l>1
        Dist2(s)=max((xim(1:l-1)-xim(2:l)).^2+(yim(1:l-1)-yim(2:l)).^2);
     end

     im_4(s,1:4)= [min(xim) max(xim) min(yim) max(yim)];

 end

Xre=cell(1,lre);
Yre=cell(1,lre);
Nre=zeros(1,lre);
Dist1=zeros(1,lim);
re_4=zeros(lim,4);
 for s= 1:lre

     xre =chre(s).x;
     yre =chre(s).y;
     l=length(xre);

     Xre{s} =xre(1:l);
     Yre{s} =yre(1:l);  
     Nre(s)     = l;

     if l>1
        Dist1(s)=max((xre(1:l-1)-xre(2:l)).^2+(yre(1:l-1)-yre(2:l)).^2);
     end

     re_4(s,1:4)= [min(xre) max(xre) min(yre) max(yre)];

 end


 %% find out which contour speaks to which
 M = zeros(lre,lim);

 for t =1:lre

     min_x_re = re_4(t,1);
     max_x_re = re_4(t,2);
     min_y_re = re_4(t,3);
     max_y_re = re_4(t,4);

     for s=1:lim

         min_x_im = im_4(s,1);
         max_x_im = im_4(s,2);
         min_y_im = im_4(s,3);
         max_y_im = im_4(s,4);

         if ( (((min_x_re <= min_x_im) && (min_x_im <= max_x_re)) || ...
               ((min_x_re <= max_x_im) && (max_x_im <= max_x_re)) || ...
               ((min_x_im <= min_x_re) && (min_x_re <= max_x_im)) || ...
               ((min_x_im <= max_x_re) && (max_x_re <= max_x_im))) && ...
              (((min_y_re <= min_y_im) && (min_y_im <= max_y_re)) || ...
               ((min_y_re <= max_y_im) && (max_y_im <= max_y_re)) || ...
               ((min_y_im <= min_y_re) && (min_y_re <= max_y_im)) || ...
               ((min_y_im <= max_y_re) && (max_y_re <= max_y_im))) )

             M(t,s)=1;
         end

     end
 end


%% locate pinwheel position
for t=1:lre

  lr=Nre(t);
  xre = Xre{t};
  yre = Yre{t};

  dist1=Dist1(t);

  mm  = M(t,:);
  s2t = find(mm==1);

  for ss= 1:length(s2t);

     s = s2t(ss);

     li   = Nim(s);

     xim = Xim{s};
     yim = Yim{s};

     dist2=Dist2(s);

     de = 2*max(dist1,dist2);

     
    [idx,~]= rangesearch([xim;yim]',[xre;yre]',sqrt(de),'NSMethod','kdtree','Distance','euclidean');
     for p=1:lr-1

%         dist = (xre(p)-xim).^2 + (yre(p)-yim).^2;
%         index= find (dist <= de);
        
        index=idx{p};
        for q=1:length(index);

            if index(q)<=li-1;

                cx=xim(index(q));
                cy=yim(index(q));

                dx=xim(index(q)+1);
                dy=yim(index(q)+1);
                ax=xre(p);
                ay=yre(p);
                bx=xre(p+1);
                by=yre(p+1);

                [answer,m,~] =  dotheycross(ax,ay,bx,by,cx,cy,dx,dy);
                if (answer)

                    pinx=cx+m*(dx-cx);
                    piny=cy+m*(dy-cy);
                    
                    x=round(pinx);
                    y=round(piny);

                    %checking whether pinwheel is at boundary
                    if x>1 && x < sx && y>1 && y <sy
                        det=grx(y,x)*giy(y,x)-gry(y,x)*gix(y,x);  % pinwheel sign (notice transpose)
                        % check for pinwheels counted twice
                        if   ~((x==doppeltestx) && (y==doppeltesty) && (sign(det)==doppeldet))
                            count=count+1;
                            PWxList=[PWxList pinx];
                            PWyList=[PWyList piny];
                            signList=[signList det];

                            doppeltestx=x;
                            doppeltesty=y;
                            doppeldet=sign(det);

                        end
                    end

                 end

           end
         end


      end


  end


 end






end


function [answer,m,n] = dotheycross(ax,ay,bx,by,cx,cy,dx,dy)

 m=((cx-ax)*(by-ay)-(cy-ay)*(bx-ax))/((cx-dx)*(by-ay)-(cy-dy)*(bx-ax));
 n=((cx-ax)*(cy-dy)-(cy-ay)*(cx-dx))/((bx-ax)*(cy-dy)-(by-ay)*(cx-dx));

 answer=((0 <= m) && (m <= 1) && (0<=n) &&(n<=1));

end

 
 function s = getcontourlines(c)

    sz = size(c,2);     % Size of the contour matrix c
    ii = 1;             % Index to keep track of current location
    jj = 1;             % Counter to keep track of # of contour lines

    while ii < sz       % While we haven't exhausted the array
        n = c(2,ii);    % How many points in this contour?
        s(jj).x = c(1,ii+1:ii+n); % X coordinates
        s(jj).y = c(2,ii+1:ii+n); % Y coordinates
        ii = ii + n + 1;          % Skip ahead to next contour line
        jj = jj + 1;              % Increment number of contours
    end

 end


