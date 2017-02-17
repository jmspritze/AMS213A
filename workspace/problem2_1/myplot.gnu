      # Gnuplot script file for plotting data in file "output.txt"
      # This file is called myplot.gnu
      set   autoscale                        # scale axes automatically
      unset log                              # remove any log-scaling
      unset label                            # remove any previous labels
      set xtic auto                          # set xtics automatically
      set ytic auto                          # set ytics automatically
      set title "Cholesky LS - Order 3 & 5"
      set xlabel "dx"
      set ylabel "b/Solutions"
     # set key 0.01,100
     # set label "Yield Point" at 0.003,260
     # set arrow from 0.0028,250 to 0.003,280
      set xr [0.0:1]
      set yr [0:2]
      plot    "output.txt" using 1:2 title 'O3 Calc points' with linespoints , \
              "output5.txt" using 1:2 title 'O5 Calc points' with linespoints , \
            "atkinson.dat" using 1:2 title 'Orig points' with points lc rgb "red"

      set terminal png
      set output "output.png"
      
