pro idljobs

readcol,'Ka1_c1.00000_var_window.txt',l,b
readcol,'../data/wmap/var60_wmap_Ka1_ampl_bl_5yr_v3.txt',l2,b2

plot,l,b
oplot,l2,b2,psym=2


  
end
