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

```
# Other randoms things that I am trying
sudo apt install gcc-8 g++-8 gcc-9 g++-9 gcc-10 g++-10
sudo cpan PerlX::Window

Running make test for TOBYINK/PerlX-Window-0.004.tar.gz
PERL_DL_NONLAZY=1 "/usr/bin/perl" "-MExtUtils::Command::MM" "-MTest::Harness" "-e" "undef *Test::Harness::Switches; test_harness(0, 'blib/lib', 'blib/arch')" t/*.t
t/01basic.t .. Data::Alias confused in da_ck_entersub (da_inside < 0) at /root/.cpan/build/PerlX-Window-0.004-1/blib/lib/PerlX/Window.pm line 135.
Compilation failed in require at t/01basic.t line 27.
BEGIN failed--compilation aborted at t/01basic.t line 27.
t/01basic.t .. Dubious, test returned 255 (wstat 65280, 0xff00)
No subtests run 
```
