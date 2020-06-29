function dbg(v)
vp = csvread('../var.csv');

dv = v(:) - vp(:);

if max(sum(abs(dv)))<1e-7
     cprintf('_green', 'Pass');
else
    
    dv = v(:) - (vp(:)+1);
    if max(sum(abs(dv)))<1e-7
         cprintf('_blue', 'Pass with +1 of python');
    else
     cprintf('_red', 'Failed');
    end
end


end