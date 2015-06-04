function [Y, FWHM, integral, X] = bellCurve(a,b,c,N)
%function [Y FWHM integral X] = bellCurve(a,b,c,N)

X=linspace(-N+b,N+b,N);
Y = a * exp( -((X-b).^2/(2*c^2)) );

FWHM=2*sqrt(2*log(2))*c;

integral=a*c*sqrt(2*pi);

plotIt=0;
if plotIt==1
    plot(X,Y,'b')
    hold on
    plot([b b],[0 a],'r-')
    plot([b-FWHM/2 b+FWHM/2],[a/2 a/2],'g-')
    plot([b-FWHM/2 b-FWHM/2],[0 a/2],'g:')
    plot([b+FWHM/2 b+FWHM/2],[0 a/2],'g:')
    plot([b-FWHM/2 b+FWHM/2],[0 0],'g-')
    hold off
    box off
end

