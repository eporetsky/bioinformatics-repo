# Gene colocality analysis

This will be an attempt to reproduce the following paper:\
https://www.osti.gov/pages/biblio/1659687\
https://github.com/ffoflonker/gene-neighborhoods

# Content

Run the following command to get all the Phytozome files from the nested folder into one folder.
```
Should work from within the base unzipped folder
cp */*/*/* ../
```

Use excel text-to-column with "ls -l" to relatively quickly get the species names from files.

Run the following sed command to batch rename all files to just species name
```
Use a conversion csv file in the following form: "from.gz,to.gz" 
Wonderful explanation here: https://unix.stackexchange.com/questions/57754/how-to-rename-files-with-sed-and-csv
Also looks like Phytozome has DOS-style line endings (\r\n) and the \r needs to be removed as well.

sed 's/^/mv /;s/,/ /;s/\r//;' < conversion.csv | bash - 
```