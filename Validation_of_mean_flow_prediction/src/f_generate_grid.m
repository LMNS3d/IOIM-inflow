
% Grid generation (adapted from Xcompact3d)

function yp = f_generate_grid(ny,beta1,yly)

    %%%
    % Alloc and init
    yeta = zeros(ny,1);

    % Local value for nym
    nym = ny-1;

    yinf = -yly/2;
    den  = 2*beta1*yinf;
    xnum = -yinf-sqrt(pi*pi*beta1*beta1+yinf*yinf);
    alpha=abs(xnum/den);

    yp = zeros(ny,1);

    if alpha ~= 0.
       yp(1)   = 0;
       yeta(1) = -0.5;
       for j=2:ny
          yeta(j)=real(j-1)*(0.5/nym)-0.5;
          den1=sqrt(alpha*beta1+1);
          xnum=den1/sqrt(alpha/pi)/sqrt(beta1)/sqrt(pi);
          den=2*sqrt(alpha/pi)*sqrt(beta1)*pi*sqrt(pi);
          den3=((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta1/pi)+alpha/pi;
          den4=2*alpha*beta1-cos(2*pi*yeta(j))+1;
          xnum1=(atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den;
          cst=sqrt(beta1)*pi/(2*sqrt(alpha)*sqrt(alpha*beta1+1));
          if yeta(j)< 0.5
              yp(j)=(xnum1-cst+yly)*2;
          end
          if yeta(j) == 0.5
              yp(j)=(0+yly)*2;
          end
          if yeta(j) > 0.5
              yp(j)=(xnum1+cst+yly)*2;
          end
        end
    end

end

