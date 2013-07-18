cat ../2.chis/*chi | awk '/^pilot.wendao.cenrot/ { print $5, $6, $7 > $3".dat"}'
