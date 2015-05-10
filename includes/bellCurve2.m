function [G, FWHM] = bellCurve2(height,centers,sigmas,imSize,theta)

[X, Y]=meshgrid(1:imSize(1),1:imSize(2));

muX=centers(1);
muY=centers(2);
sigmaX=sigmas(1);
sigmaY=sigmas(2);

a=cos(theta)^2 /(2*sigmaX^2) + sin(theta)^2 /(2*sigmaY^2);
b=-sin(2*theta) /(4*sigmaX^2) + sin(2*theta) /(4*sigmaY^2);
c=sin(theta)^2 /(2*sigmaX^2) + cos(theta)^2 /(2*sigmaY^2);

G=height * exp( -(a*(X-muX).^2 + 2*b*(X-muX).*(Y-muY) + c*(Y-muY).^2 ) );

[x, FWHM1]=bellCurve(height,muX,sigmaX,imSize(1));
[x, FWHM2]=bellCurve(height,muX,sigmaX,imSize(1));
FWHM=[FWHM1 FWHM2];

plotIt=0;
if plotIt==1
    imshow(G,[])
    title(sci([[a b c]*1000 theta]))
end


