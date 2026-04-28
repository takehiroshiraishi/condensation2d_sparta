set terminal x11 0 persist
set xlabel "column 1"
set ylabel "column 3"
set key outside
set xrange [0:*]

title(d) = system(sprintf("echo '%s' | awk -F'hbox_' '{v=$2; gsub(/p/,\".\",v); printf \"%.0fum\", v*1e6}'", d))

plot for [d in system("awk -F/ '{print $NF}' case_list.txt")] \
    sprintf("%s/profiles_steady/y_axis.dat", d) u 1:3 w l title title(d)


set terminal x11 1 persist
set xlabel "column 1"
set ylabel "column 2"
set key outside
set xrange [0:*]

plot for [d in system("awk -F/ '{print $NF}' case_list.txt")] \
    sprintf("%s/profiles_steady/y_axis.dat", d) u 1:2 w l title title(d)
