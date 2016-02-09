do for [i=0:40] {
  plot "force.out" every ::i::(i+1)*1600 using 1:2:3 with points palette pointsize 3 pointtype 7 title "w='i'"
}
