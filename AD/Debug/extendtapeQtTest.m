for i=100000:110000
   if(sum(not(TAPE_old_way(i,:) == TAPE_old_way(i,:))))
       i
   end    
end