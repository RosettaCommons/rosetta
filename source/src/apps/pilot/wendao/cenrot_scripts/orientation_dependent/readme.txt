run general_counting to get load_table m file
load angle-r data and generate orientation dependent term

cat score_range.dat | awk '{print $2}' | sed 's/=/ /' | awk '{print $2}' | sort -n > dat.max
cat score_range.dat | awk '{print $3}' | sed 's/=/ /' | awk '{print $2}' | sort -nr > dat.min
