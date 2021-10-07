cat md5.csv | 
  while read line 
  do 
    FILENAME=$(echo $line | cut -d',' -f1 )
    REALMD5SUM=$(echo $line | cut -d',' -f2 )
    CALCMD5SUM=$(md5sum $FILENAME | cut -d' ' -f1)
    if [ $REALMD5SUM = $CALCMD5SUM ]; then 
      echo $FILENAME " is correct"   
    else
      echo $FILENAME " is incorrect"
    fi
done

