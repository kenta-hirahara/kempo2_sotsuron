%***************************************************************
% fft of single real function
%   input data have 2n elements
%     isign = 1 for fourier transform
%          transformd data should be multiplied by 1/n
%          with sequance of
%          c(0),c(n),c(1),s(1),c(2),s(2),......c(n-1),s(n-1)
%       while x(j) : j=1,2,....,2n is expressed as
%          x(j)= 0.5*c(0) + sum(c(k)cos(..)+s(k)sin(..)) + 0.5*c(n)
%
%     isign = -1 for inverse fourier transform
%     
%     reference: numerical recipes by w.h. press et al., cambridge 1986
%        modified by y. omura, september, 1989
%
      function[data] = realft(data,n2,isign)
      n=n2/2;
      phase=pi/n;
      c1=0.5;
      if isign == 1 
        c2=-0.5;
        data = four1(data,n2,1);
      else
        c2=0.5;
        phase=-phase;
      end
      wpr=-2.d0*sin(0.5*phase)^2;
      wpi=sin(phase);
      wr=1.0+wpr;
      wi=wpi;
      n2p3=2*n+3;
      for i1 = 3:2:n-1 
        wrs=wr;
        wis=wi;
        h1r=c1*(data(i1)+data(n2p3-i1-1));
        h1i=c1*(data(i1+1)-data(n2p3-i1));
        h2r=-c2*(data(i1+1)+data(n2p3-i1));
        h2i=c2*(data(i1)-data(n2p3-i1-1));
        data(i1)=h1r+wrs*h2r-wis*h2i;
        data(i1+1)=h1i+wrs*h2i+wis*h2r;
        data(n2p3-i1-1)=h1r-wrs*h2r+wis*h2i;
        data(n2p3-i1)=-h1i+wrs*h2i+wis*h2r;
        wtemp=wr;
        wr=wr*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
      end
      if isign == 1 
        h1r=data(1);
        data(1)=h1r+data(2);
        data(2)=h1r-data(2);
      else
        h1r=data(1);
        data(1)=c1*(h1r+data(2));
        data(2)=c1*(h1r-data(2));
        data = four1(data,n2,-1);
      end
  return
    

