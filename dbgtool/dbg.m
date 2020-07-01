function dbg(v)
clc
if iscell(v)
   for i =1:length(v)
      vi = v{i};
      if iscomplex(vi)          
          fn = sprintf('C:/Users/yansh/Desktop/CanonicalCoordinateGeometryLearn/%d_real.csv',i);
          basic_check(fn,real(vi))    
          fn = sprintf('C:/Users/yansh/Desktop/CanonicalCoordinateGeometryLearn/%d_imag.csv',i);
          basic_check(fn,imag(vi))    
      else
          fn = sprintf('C:/Users/yansh/Desktop/CanonicalCoordinateGeometryLearn/%d.csv',i);
          basic_check(fn,vi)    
      end
      
       
   end
    
else
    
    if iscomplex(v)
        fn = 'C:/Users/yansh/Desktop/CanonicalCoordinateGeometryLearn/real.csv';
        basic_check(fn,real(v))
        fn = 'C:/Users/yansh/Desktop/CanonicalCoordinateGeometryLearn/imag.csv';
        basic_check(fn,imag(v))
    else
        fn = 'C:/Users/yansh/Desktop/CanonicalCoordinateGeometryLearn/var.csv';
        basic_check(fn,v)
    end
end



end


function flag = basic_check(fn, v)

vp = csvread(fn);
dv = v(:) - vp(:);
flag = -1;
if max(abs(dv))<1e-7
     cprintf('_green', 'Pass');
     flag = 0;
else
    
    dv = v(:) - (vp(:)+1);
    if max(abs(dv))<1e-7
         cprintf('_blue', 'Pass with +1 of python');
         flag=1;
    else
        
     [max1, i1]=max(abs(v-vp));
     [max2, i2]=max(max1);
     
     
        
     cprintf('_red', 'Failed');
     
     cprintf('_red', sprintf('i=%d, j = %d with maxdiff=%10.10f/%10.10f\n', i1(i2),i2, max2, v(i1(i2))));
     flag = -1;
    end
end

end