set xlabel 'eps_v'
set ylabel 'sig_v'
sv(sxx,syy,szz,sxy,syz,sxz)=sqrt(  0.5*( (sxx-syy)**2 + (syy-szz)**2 + (sxx-szz)**2) + 6.*(syz**2+sxz**2+sxy**2)  )


# sv=sqrt(  0.5*( (sxx-syy).^2+(syy-szz).^2+(sxx-szz).^2) + 6.*(syz.^2+sxz.^2+sxy.^2));
p 'task2.in.out.hom' u (sv($2,$3,$4,$5,$6,$7)) : (sv($8,$9,$10,$11,$12,$13))  w l t'elastic'

pause -1

p 'task2.in.out.hom' u (sv($2,$3,$4,$5,$6,$7)) *(sv($8,$9,$10,$11,$12,$13))  w l t'total (av.) strain energy density'
pause -1
