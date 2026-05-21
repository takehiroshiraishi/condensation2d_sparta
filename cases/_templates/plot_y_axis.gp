set terminal x11 0 persist
set xlabel "distance [m]"
set ylabel "temperature [K]"
set key outside
set xrange [0:*]

title(d) = d

plot for [d in system("awk -F/ '{print $NF}' case_list.txt")] \
    sprintf("%s/profiles_steady/y_axis.dat", d) u 1:3 w l title title(d)


set terminal x11 1 persist
set xlabel "distance [m]"
set ylabel "pressure [Pa]"
set key outside
set xrange [0:*]

plot for [d in system("awk -F/ '{print $NF}' case_list.txt")] \
    sprintf("%s/profiles_steady/y_axis.dat", d) u 1:2 w l title title(d)


set terminal x11 2 persist
set xlabel "distance [m]"
set ylabel "mass flux in -y [kg/m^2/s]"
set key outside
set xrange [0:*]

plot for [d in system("awk -F/ '{print $NF}' case_list.txt")] \
    sprintf("%s/profiles_steady/y_axis.dat", d) u 1:5 w l title title(d)
