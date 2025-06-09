close all, clear all, clc

result = -5;

switch(result)
   case 52
      disp('result is 52')
   case {52, 78}
      disp('result is 52 or 78')
    otherwise
        disp('Donno')
end


