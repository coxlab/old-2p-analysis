function timeStr = sec2time(x)
%function timeStr = sec2time(x)

hour = floor(x/60^2);
minutes = floor(mod((x/60), 60));
seconds = floor(mod(x,60));
%miliseconds = mod(x,1)*1e3;

timeStr = sprintf('%02d:%02d:%02d', hour, minutes, seconds);
