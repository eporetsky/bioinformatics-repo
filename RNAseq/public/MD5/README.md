# Compare MD5sum between origin and downloaded file

An md5sum csv files is used to compare original vs. downloaded md5sum results\
The provided csv file (check_md5sum_test.csv) can be used to test the script.\ 
Script (check_md5sum.sh) and files should be in the same folder.\

## To run the check_md5sum.bash:
```
bash check_md5sum.sh -f check_md5sum_test.csv
```

## The bash script will iterate over the rows of the csv file
```
# The test csv file will test the md5sum of the bash file
check_md5sum.sh  is correct
```
